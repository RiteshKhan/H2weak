//*****************************************************************
// Iterative template routine -- GMRES
//
// Remodified from : https://github.com/amiraa127/Sparse_MultiFrontal/blob/master/IML/include/gmres.h
// GMRES solves the unsymmetric linear system Ax = b using the
// Generalized Minimum Residual method (GMRES)

#ifndef IML_GMRES
#define IML_GMRES
#include "kernel.hpp"
#include "matvec.hpp"
using namespace std;

class iterative_solver
{
public:
    int restart, max_iter;
    string fname;
    double gmres_tol;
    double resid;

    iterative_solver(int max_iter, int gmres_tol_pow)
    {
        this->restart = max_iter;
        this->max_iter = max_iter;
        this->gmres_tol = pow(10, -1.0 * gmres_tol_pow);
        this->resid = 0.0;
    }

    void GeneratePlaneRotation(dtype &dx, dtype &dy, dtype &cs, dtype &sn)
    {
#ifdef USE_real
        if (dy == 0.0)
        {
            cs = 1.0;
            sn = 0.0;
        }
        else if (abs(dy) > abs(dx))
        {
            double temp = dx / dy;
            sn = 1.0 / sqrt(1.0 + temp * temp);
            cs = temp * sn;
        }
        else
        {
            double temp = dy / dx;
            cs = 1.0 / sqrt(1.0 + temp * temp);
            sn = temp * cs;
        }
#endif

#ifdef USE_Helm
        if (dx == dtype(0))
        {
            cs = dtype(0);
            sn = dtype(1);
        }
        else
        {
            double scale = std::abs(dx) + std::abs(dy);
            double norm = scale * std::sqrt(std::abs(dx / scale) * std::abs(dx / scale) +
                                            std::abs(dy / scale) * std::abs(dy / scale));
            dtype alpha = dx / std::abs(dx);
            cs = std::abs(dx) / norm;
            sn = alpha * std::conj(dy) / norm;
        }
#endif
    }

    void ApplyPlaneRotation(dtype &dx, dtype &dy, dtype &cs, dtype &sn)
    {
        dtype temp = cs * dx + sn * dy;
        dy = -sn * dx + cs * dy;
        dx = temp;
    }

    void Update(Vec &x, int k, Mat &h, Vec &s, Vec *&v)
    {
        Vec y = s;
        // Backsolve:
        for (int i = k; i >= 0; i--)
        {
            y(i) /= h(i, i);
            for (int j = i - 1; j >= 0; j--)
                y(j) -= h(j, i) * y(i);
        }

        for (int j = 0; j <= k; j++)
            x += v[j] * y(j);
    }
    template <typename T>
    void GMRES(T *&A, Vec &x, Vec &b)
    {
        std::cout << "GMRES Parameters " << std::endl;
        std::cout << "Maximum Iterations : " << max_iter << std::endl;
        std::cout << "GMRES Tolerance : " << gmres_tol << std::endl;

        int m = restart;
        Mat H = Mat::Zero(restart + 1, restart);
        int i, j = 1, k;
        double gmres_mvp = 0.0;

        // Vectors for Rotation

        Vec s = Vec::Zero(restart + 1), cs = Vec::Zero(restart + 1), sn = Vec::Zero(restart + 1), w;
        double normb = b.norm();
        Vec output1;

        clock_t start_gmres_mvp1;
        start_gmres_mvp1 = clock();
        A->mat_vec_prod(x, output1);
        gmres_mvp += (clock() - start_gmres_mvp1) / (double)CLOCKS_PER_SEC;

        Vec r = b - output1; // Explict method to perform the MatVec

        double beta = r.norm();

        if (normb == 0.0)
            normb = 1;
        resid = beta / normb;
        if (resid <= gmres_tol)
        {
            gmres_tol = resid;
            max_iter = 0;
            std::cout << "Early Exited" << std::endl;
            return;
        }
        Vec *v = new Vec[restart + 1];
        while (j <= max_iter)
        {

            v[0] = r * (1.0 / beta);
            s = Vec::Zero(restart + 1);
            s(0) = beta;

            for (i = 0; i < m && j <= max_iter; i++, j++)
            {
                // std::cout << "Iteration[" << j << "] : " << std::setprecision(10) << resid << std::endl;
                Vec output2;

                clock_t start_gmres_mvp2;
                start_gmres_mvp2 = clock();
                A->mat_vec_prod(v[i], output2);
                gmres_mvp += (clock() - start_gmres_mvp2) / (double)CLOCKS_PER_SEC;

                w = output2;

                for (k = 0; k <= i; k++)
                {
                    H(k, i) = w.dot(v[k]);
                    w -= H(k, i) * v[k];
                }
                H(i + 1, i) = w.norm();
                v[i + 1] = (1.0 / H(i + 1, i)) * w;

                for (k = 0; k < i; k++)
                    ApplyPlaneRotation(H(k, i), H(k + 1, i), cs(k), sn(k));

                GeneratePlaneRotation(H(i, i), H(i + 1, i), cs(i), sn(i));
                ApplyPlaneRotation(H(i, i), H(i + 1, i), cs(i), sn(i));
                ApplyPlaneRotation(s(i), s(i + 1), cs(i), sn(i));
                resid = abs(s(i + 1)) / normb;
                if (resid < gmres_tol)
                {
                    Update(x, i, H, s, v);
                    gmres_tol = resid;
                    max_iter = j;
                    std::cout << std::endl
                              << "Reached Solution before Max_Iterations " << std::endl;
                    std::cout << "Reached the desired tol after " << j << " iterations" << std::endl;
                    std::cout << "Final resid " << resid << std::endl;
                    //std::cout << "The GMRES mvp time " << gmres_mvp << std::endl;
                    delete[] v;
                    return;
                }
            }
            std::cout << "Restarted" << std::endl;
            Update(x, m - 1, H, s, v);

            Vec output3;
            clock_t start_gmres_mvp3;
            start_gmres_mvp3 = clock();
            A->mat_vec_prod(x, output3);
            gmres_mvp += (clock() - start_gmres_mvp3) / (double)CLOCKS_PER_SEC;

            r = b - output3;

            beta = r.norm();
            resid = beta / normb;
            if (resid < gmres_tol)
            {
                max_iter = j;
                std::cout << "Reached desired tolarence before Max_Iterations" << std::endl;
                delete[] v;
                return;
            }
        }
        delete[] v;
        return;
    }
};

#endif