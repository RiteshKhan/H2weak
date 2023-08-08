#ifndef _KERNEL_HPP__
#define _KERNEL_HPP__
#include "integral2d.hpp"

class kernel
{
public:
	double a;
	std::vector<pts2D> particles;

	kernel(std::vector<pts2D> &particles)
	{
		this->particles = particles;
	}

	virtual dtype getMatrixEntry(const unsigned i, const unsigned j)
	{
		std::cout << "virtual getInteraction" << std::endl;
		return 0.0;
	}

	Vec getRow(const int j, std::vector<int> col_indices)
	{
		int n_cols = col_indices.size();
		Vec row(n_cols);
#pragma omp parallel for
		for (int k = 0; k < n_cols; k++)
		{
			row(k) = this->getMatrixEntry(j, col_indices[k]);
		}
		return row;
	}

	Vec getCol(const int k, std::vector<int> row_indices)
	{
		int n_rows = row_indices.size();
		Vec col(n_rows);
#pragma omp parallel for
		for (int j = 0; j < n_rows; ++j)
		{
			col(j) = this->getMatrixEntry(row_indices[j], k);
		}
		return col;
	}

	Vec getCol(const int k)
	{
		int n_rows = particles.size();
		Vec col(n_rows);
#pragma omp parallel for
		for (int j = 0; j < n_rows; ++j)
		{
			col(j) = this->getMatrixEntry(j, k);
		}
		return col;
	}

	Mat getMatrix(std::vector<int> row_indices, std::vector<int> col_indices)
	{
		int n_rows = row_indices.size();
		int n_cols = col_indices.size();
		Mat mat(n_rows, n_cols);
		for (int j = 0; j < n_rows; ++j)
		{
			for (int k = 0; k < n_cols; ++k)
			{
				mat(j, k) = this->getMatrixEntry(row_indices[j], col_indices[k]);
			}
		}
		return mat;
	}

	Mat getMatrix(int row_start_index, int col_start_index, int row_end_index, int col_end_index)
	{
		Mat mat(row_end_index - row_start_index, col_end_index - col_start_index);
		for (int j = row_start_index; j < row_end_index; ++j)
		{
			for (int k = col_start_index; k < col_end_index; ++k)
			{
				mat(j, k) = this->getMatrixEntry(j, k);
			}
		}
		return mat;
	}

	~kernel(){};
};

class userkernel : public kernel
{
public:
	int ker_choice;
	dtype Kii;
	double h2;
	double chargesFunction(const pts2D r)
	{
		double q = r.x; // user defined
		return q;
	};

	userkernel(std::vector<pts2D> particles, int ker_choice) : kernel(particles)
	{
		this->ker_choice = ker_choice;
		double h = 1.0 / sqrt(double(particles.size()));
		double *a, *b;
		a = new double[2];
		b = new double[2];
		a[0] = 0;
		a[1] = 0;

		b[0] = h * 0.5;
		b[1] = h * 0.5;

		Kii = double_integral(a, b);
		h2 = pow(h, 2.0);
	};

	dtype IE_2D(const unsigned i, const unsigned j)
	{
		if (i == j)
			return (1.0 + Kii) / (double) h2; // For Second kind integral equation.
		else
			return Laplacian_2D(i, j);
	}

#ifdef USE_Hankel

	dtype IE_2D_Helm(const unsigned i, const unsigned j)
	{
		if (i == j)
			return (1.0 + 15*log(double(particles.size())) + Kii) / (double) h2 ; // For Second kind integral equation.
		else
			return Helmholtz2D(i, j);
	}

	dtype Helmholtz2D(const unsigned i, const unsigned j)
	{
		pts2D r1 = particles[i]; // particles_X is a member of base class FMM_Matrix
		pts2D r2 = particles[j]; // particles_X is a member of base class FMM_Matrix
		double R2 = (r1.x - r2.x) * (r1.x - r2.x) + (r1.y - r2.y) * (r1.y - r2.y);
		double R = sqrt(R2);

		return I * (besselJ(0, R) + I * besselY(0, R)) / 4.0;
	}
#endif

	double RBF_Logarithm(const unsigned i, const unsigned j)
	{
		pts2D r1 = particles[i];
		pts2D r2 = particles[j];
		double R2 = (r1.x - r2.x) * (r1.x - r2.x) + (r1.y - r2.y) * (r1.y - r2.y);
		double R = sqrt(R2);
		double b = 1.0;
		return log(1.0 + R / b);
	}

	double RBF_Exponential(const unsigned i, const unsigned j)
	{
		pts2D r1 = particles[i];
		pts2D r2 = particles[j];
		double R2 = (r1.x - r2.x) * (r1.x - r2.x) + (r1.y - r2.y) * (r1.y - r2.y);
		double R = sqrt(R2);
		double b = 1.0;
		return exp(-R / b);
	}

	double RBF_Inverse_Quadratic(const unsigned i, const unsigned j)
	{
		pts2D r1 = particles[i];
		pts2D r2 = particles[j];
		double R2 = (r1.x - r2.x) * (r1.x - r2.x) + (r1.y - r2.y) * (r1.y - r2.y);
		double R = sqrt(R2);
		double b = 1.0;
		return 1.0 / (1.0 + (R / b) * (R / b));
	}

	double RBF_Sin(const unsigned i, const unsigned j)
	{
		pts2D r1 = particles[i];
		pts2D r2 = particles[j];
		double R2 = (r1.x - r2.x) * (r1.x - r2.x) + (r1.y - r2.y) * (r1.y - r2.y);
		double R = sqrt(R2);
		double b = 1.0;
		return sin(R / b);
	}

	double RBF_Gaussian(const unsigned i, const unsigned j)
	{
		pts2D r1 = particles[i];
		pts2D r2 = particles[j];
		double R2 = (r1.x - r2.x) * (r1.x - r2.x) + (r1.y - r2.y) * (r1.y - r2.y);
		double R = sqrt(R2);
		double b = 1.0;
		return exp(-(R / b) * (R / b));
	}

	double R_over_A_over_R(const unsigned i, const unsigned j)
	{
		pts2D r1 = particles[i];
		pts2D r2 = particles[j];
		double R2 = (r1.x - r2.x) * (r1.x - r2.x) + (r1.y - r2.y) * (r1.y - r2.y);
		double R = sqrt(R2);
		double a = 0.0001;
		if (R < a)
		{
			return R / a;
		}
		else
		{
			return a / R;
		}
	}

	double Laplacian_2D(const unsigned i, const unsigned j)
	{
		pts2D r1 = particles[i];
		pts2D r2 = particles[j];
		double R2 = (r1.x - r2.x) * (r1.x - r2.x) + (r1.y - r2.y) * (r1.y - r2.y);
		double c = -1.0 / (2 * PI);
		if (i == j)
		{
			return 0.0;
		}
		else
		{
			return c * 0.5 * log(R2);
		}
	}

	double Laplacian_3D(const unsigned i, const unsigned j)
	{
		pts2D r1 = particles[i];
		pts2D r2 = particles[j];
		double R2 = (r1.x - r2.x) * (r1.x - r2.x) + (r1.y - r2.y) * (r1.y - r2.y);
		double R = sqrt(R2);
		if (i == j)
		{
			return 0.0;
		}
		else
		{
			return 1.0 / R;
		}
	}

	// 1/r^4
	double oneOverR4(const unsigned i, const unsigned j)
	{
		pts2D r1 = particles[i];
		pts2D r2 = particles[j];
		double R2 = (r1.x - r2.x) * (r1.x - r2.x) + (r1.y - r2.y) * (r1.y - r2.y);
		double R = sqrt(R2);
		double b = 0.001;
		if (R < b)
		{
			return pow(R / b, 4);
		}
		else
		{
			return pow(b / R, 4);
		}
	}

	// RBF r^{2}log(r)
	double RBF_spline(const unsigned i, const unsigned j)
	{
		pts2D r1 = particles[i];
		pts2D r2 = particles[j];
		double R2 = (r1.x - r2.x) * (r1.x - r2.x) + (r1.y - r2.y) * (r1.y - r2.y);
		double R = sqrt(R2);
		if (i == j)
		{
			return 0;
		}
		else
		{
			return R * R * log(R);
		}
	}

	double Helmholtz_cos(const unsigned i, const unsigned j)
	{
		pts2D r1 = particles[i];
		pts2D r2 = particles[j];
		double R2 = (r1.x - r2.x) * (r1.x - r2.x) + (r1.y - r2.y) * (r1.y - r2.y);
		double R = sqrt(R2);
		// double a = 0.001;
		double kappa = 1.0;
		if (i == j)
		{
			return 0;
		}
		else
		{
			return cos(kappa * R) / R;
		}
	}
#ifdef USE_real
	dtype getMatrixEntry(const unsigned i, const unsigned j)
	{
		if (ker_choice != 0)
		{
			std::cout << "Wrong choice! Please give correct input " << std::endl;
			exit(1);
		}
		else
		{
			if (ker_choice == 0)
				return IE_2D(i, j);
			if (ker_choice == 1)
				return RBF_Logarithm(i, j);
			else if (ker_choice == 2)
				return RBF_Exponential(i, j);
			else if (ker_choice == 3)
				return RBF_Inverse_Quadratic(i, j);
			else if (ker_choice == 4)
				return RBF_Sin(i, j);
			else if (ker_choice == 5)
				return RBF_Gaussian(i, j);
			else if (ker_choice == 6)
				return R_over_A_over_R(i, j);
			else if (ker_choice == 7)
			{
				if (i == j)
					return 10000.0;
				else
					return Laplacian_2D(i, j);
			}
			else if (ker_choice == 8)
			{
				if (i == j)
					return 0.0;
				else
					return Laplacian_3D(i, j);
			}
			else if (ker_choice == 9)
			{
				if (i == j)
					return 0.0;
				else
					return oneOverR4(i, j);
			}
			else if (ker_choice == 10)
			{
				return RBF_spline(i, j);
			}
			else
			{
				return Helmholtz_cos(i, j);
			}
		}
	}
#endif

#ifdef USE_Hankel
	dtype getMatrixEntry(const unsigned i, const unsigned j)
	{
		if (ker_choice == 10)
		{
			return IE_2D_Helm(i, j);
		}
		else
		{
			std::cout << "Wrong choice! Please give correct input " << std::endl;
			exit(1);
		}
	}
#endif

	~userkernel(){};
};

#endif