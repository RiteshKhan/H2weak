#ifndef _KERNEL_HPP__
#define _KERNEL_HPP__
#include "integral3d.hpp"

class kernel
{
public:
	double a;
	std::vector<pts3D> particles;

	kernel(std::vector<pts3D> &particles)
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
	int ker_choice, N;
	dtype Kii;
	double h3;
	double chargesFunction(const pts3D r)
	{
		double q = r.x; // user defined
		return q;
	};

	userkernel(std::vector<pts3D> particles, int ker_choice) : kernel(particles)
	{
		this->ker_choice = ker_choice;
		if (ker_choice == 0)
		{
			double h = 1.0 / cbrt(double(particles.size()));
			double *a, *b;
			a = new double[3];
			b = new double[3];
			a[0] = 0;
			a[1] = 0;
			a[2] = 0;

			b[0] = h * 0.5;
			b[1] = h * 0.5;
			b[2] = h * 0.5;

			Kii = triple_integral(a, b);
			h3 = pow(h, 3.0);
		}
	};
#ifdef USE_real
	dtype IE_3D(const unsigned i, const unsigned j)
	{
		if (i == j)
			return (1.0 + Kii) / (double)h3; // For Second kind integral equation.
		else
			return Laplacian_3D(i, j);
	}
#endif

#ifdef USE_Helm
	dtype Helmholtz3D(const unsigned i, const unsigned j)
	{
		pts3D r1 = particles[i]; 
		pts3D r2 = particles[j];
		double R2 = (r1.x - r2.x) * (r1.x - r2.x) + (r1.y - r2.y) * (r1.y - r2.y) + (r1.z - r2.z) * (r1.z - r2.z);
		double R = sqrt(R2);
		if (i == j)
			return 500 * sqrt(double(particles.size())); // For Second kind integral equation.
		else
			return (exp(I * R)) / (4.0 * PI);
	}
#endif

	double LOGR_over_LOGA_over_LOGR(const unsigned i, const unsigned j)
	{
		pts3D r1 = particles[i];
		pts3D r2 = particles[j];
		double R2 = (r1.x - r2.x) * (r1.x - r2.x) + (r1.y - r2.y) * (r1.y - r2.y) + (r1.z - r2.z) * (r1.z - r2.z);
		double R = sqrt(R2);
		double a = 0.0001;
		if (R < a)
		{
			return (R * log(R) - 1.0) / (a * log(a) - 1.0);
		}
		else
		{
			return (1.0 * log(R)) / log(a);
		}
	}

	double Laplacian_3D(const unsigned i, const unsigned j)
	{
		pts3D r1 = particles[i];
		pts3D r2 = particles[j];
		double R2 = (r1.x - r2.x) * (r1.x - r2.x) + (r1.y - r2.y) * (r1.y - r2.y) + (r1.z - r2.z) * (r1.z - r2.z);
		double R = sqrt(R2);
		return 1.0 / (4.0 * PI * R);
	}

#ifdef USE_real
	dtype getMatrixEntry(const unsigned i, const unsigned j)
	{
		if (ker_choice != 0 && ker_choice != 1)
		{
			std::cout << "Wrong choice! Please give correct input " << std::endl;
			exit(1);
		}
		else
		{
			if (ker_choice == 0)
				return IE_3D(i, j);
			else
			{
				if (i == j)
					return (double) sqrt(particles.size());
				else
					return LOGR_over_LOGA_over_LOGR(i, j);
			}
		}
	}
#endif

#ifdef USE_Helm
	dtype getMatrixEntry(const unsigned i, const unsigned j)
	{
		if (ker_choice == 2)
		{
			return Helmholtz3D(i, j);
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
