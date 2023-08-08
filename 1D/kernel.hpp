#ifndef _KERNEL_HPP__
#define _KERNEL_HPP__
#include "integral1d.hpp"

class kernel
{
public:
	double a;
	std::vector<double> particles;

	kernel(std::vector<double> &particles)
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
	double h;
	userkernel(std::vector<double> particles, int ker_choice) : kernel(particles){
		this->ker_choice = ker_choice;
		h = 1.0 / double(particles.size());
		double *a, *b;
		a = new double[1];
		b = new double[1];
		a[0] = 0;

		b[0] = h * 0.5;

		Kii = single_integral(a, b);
	};

	dtype IE_1D(const unsigned i, const unsigned j)
	{
		if (i == j)
			return (1.0 + Kii) / (double) h; // For Second kind integral equation.
		else
			return Log_R(i, j);
	}

	double RBF_Exponential(const unsigned i, const unsigned j)
	{
		double r1 = particles[i];
		double r2 = particles[j];
		double R = fabs(r1-r2);
		double b = 1.0;
		return exp(-R / b);
	}

	double RBF_Inverse_Quadratic(const unsigned i, const unsigned j)
	{
		double r1 = particles[i];
		double r2 = particles[j];
		double R = fabs(r1-r2);
		double b = 1.0;
		return 1.0 / (1.0 + (R / b) * (R / b));
	}


	double R_over_A_over_R(const unsigned i, const unsigned j)
	{
		double r1 = particles[i];
		double r2 = particles[j];
		double R = fabs(r1-r2);
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

	double Log_R(const unsigned i, const unsigned j)
	{
		double r1 = particles[i];
		double r2 = particles[j];
		double R = fabs(r1-r2);
		if (i == j)
		{
			return 0.0;
		}
		else
		{
			return log(R);
		}
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
				return IE_1D(i, j);
			else{
				if(i==j)
					return sqrt(particles.size());
				else
					return R_over_A_over_R(i, j);
			}
		}
	}
#endif

	// #endif
	~userkernel(){};
};
#endif