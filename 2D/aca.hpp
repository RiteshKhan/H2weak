//
//  This ACA routine was created by Vaishnavi Gujjula
//
#ifndef __ACA_HPP__
#define __ACA_HPP__
#include "kernel.hpp"
class LowRank
{
public:
	userkernel *K;
	double tol_ACA;
	std::vector<int> row_indices;
	std::vector<int> col_indices;

	LowRank(userkernel *K, int tol_pow, std::vector<int> &row_indices, std::vector<int> &col_indices)
	{
		this->K = K;
		this->tol_ACA = pow(10, -1.0 * tol_pow);
		this->row_indices = row_indices;
		this->col_indices = col_indices;
	}

	void maxAbsVector(const Vec &v, const std::set<int> &allowed_indices,
					  dtype &max, int &index)
	{
		std::set<int>::iterator it;
		index = *allowed_indices.begin();
		max = abs(v(index));

		for (it = allowed_indices.begin(); it != allowed_indices.end(); it++)
		{
			if (abs(v(*it)) > abs(max))
			{
				index = *it;
				max = abs(v(index));
			}
		}
	}

	// void ACA_only_nodesUV(std::vector<int> &row_bases, std::vector<int> &col_bases, int &computed_rank, Mat &U, Mat &Vt)
	// {
	// 	int col_index;
	// 	int row_index;
	// 	int N1 = row_indices.size();
	// 	int N2 = col_indices.size();
	// 	Vec row(N2), col(N1), v(N2), u(N1);
	// 	Vec row_temp, col_temp;
	// 	// std::vector<Vec> AcVec;
	// 	// std::vector<Vec> ArVec;

	// 	std::vector<Vec> Uvec;
	// 	std::vector<Vec> Vvec;
	// 	computed_rank = 0;
	// 	dtype max;
	// 	int min = N1;
	// 	if (N1 > N2)
	// 	{
	// 		min = N2;
	// 	}

	// 	std::set<int> remaining_row_ind;
	// 	std::set<int> remaining_col_ind;

	// 	for (int l = 0; l < N1; l++)
	// 	{
	// 		remaining_row_ind.insert(l);
	// 	}
	// 	for (int l = 0; l < N2; l++)
	// 	{
	// 		remaining_col_ind.insert(l);
	// 	}
	// 	int l_local;
	// 	for (l_local = 0; l_local < N1; l_local++)
	// 	{
	// 		row_index = l_local;
	// 		row = K->getRow(row_indices[row_index], col_indices);
	// 		this->maxAbsVector(row, remaining_col_ind, max, col_index);
	// 		if (abs(row(int(col_index))) > pow(10, -36))
	// 		{
	// 			break;
	// 		}
	// 	}
	// 	if (l_local == N1)
	// 	{
	// 		return;
	// 	}
	// 	v = row;
	// 	row_bases.push_back(row_index);
	// 	col_bases.push_back(int(col_index));
	// 	col = K->getCol(col_indices[col_index], row_indices);
	// 	u = col / row(int(col_index));
	// 	Uvec.push_back(u);
	// 	Vvec.push_back(v);
	// 	remaining_col_ind.erase(col_index);
	// 	remaining_row_ind.erase(row_index);
	// 	computed_rank = 1;

	// 	double normS = 0.0;
	// 	double prev_normS = __DBL_MAX__;
	// 	this->maxAbsVector(col, remaining_row_ind, max, row_index);

	// 	while (u.norm() * v.norm() > tol_ACA * normS && computed_rank < min)
	// 	{
	// 		row_bases.push_back(int(row_index));
	// 		row = K->getRow(row_indices[row_index], col_indices);
	// 		row_temp = row;
	// 		for (int l = 0; l <= computed_rank - 1; ++l)
	// 		{
	// 			for (auto s = 0; s < Vvec[l].size(); s++)
	// 			{
	// 				row(s) -= Uvec[l](row_index) * Vvec[l](s);
	// 			}
	// 		}
	// 		this->maxAbsVector(row, remaining_col_ind, max, col_index);
	// 		col_bases.push_back(int(col_index));
	// 		v = row;
	// 		col = K->getCol(col_indices[col_index], row_indices);
	// 		col_temp = col;
	// 		for (int l = 0; l <= computed_rank - 1; ++l)
	// 		{
	// 			for (auto s = 0; s < Uvec[l].size(); s++)
	// 			{
	// 				col(s) -= Vvec[l](col_index) * Uvec[l](s);
	// 			}
	// 		}
	// 		u = col / row(int(col_index));
	// 		if (u.norm() < 1e-36 || v.norm() < 1e-36)
	// 		{
	// 			row_bases.pop_back();
	// 			col_bases.pop_back();
	// 			break;
	// 		}
	// 		Uvec.push_back(u);
	// 		Vvec.push_back(v);

	// 		++computed_rank;
	// 		remaining_col_ind.erase(col_index);
	// 		remaining_row_ind.erase(row_index);
	// 		if (computed_rank != 2)
	// 			prev_normS = normS;
	// 		normS = pow(normS, 2) + pow(u.norm() * v.norm(), 2);
	// 		for (int l = 0; l <= computed_rank - 1; ++l)
	// 		{
	// 			normS += 2 * abs(Uvec[l].dot(u) * Vvec[l].dot(v));
	// 		}
	// 		normS = sqrt(normS);
	// 		this->maxAbsVector(col, remaining_row_ind, max, row_index);
	// 	}
	// 	U = Mat(N1, computed_rank);
	// 	Vt = Mat(computed_rank, N2);
	// 	for (auto i = 0; i < computed_rank; i++)
	// 	{
	// 		U.col(i) = Uvec[i];
	// 		Vt.row(i) = Vvec[i];
	// 	}
	// }

	void ACA_only_nodesCUR(std::vector<int> &row_bases, std::vector<int> &col_bases, int &computed_rank, Mat &L, Mat &R)
	{
		int col_index;
		int row_index;
		int N1 = row_indices.size();
		int N2 = col_indices.size();
		Vec row(N2), col(N1), v(N2), u(N1);
		Vec row_temp, col_temp;

		std::vector<Vec> Uvec;
		std::vector<Vec> Vvec;
		computed_rank = 0;
		dtype max;
		int min = N1;
		if (N1 > N2)
		{
			min = N2;
		}

		std::set<int> remaining_row_ind;
		std::set<int> remaining_col_ind;

		for (int l = 0; l < N1; l++)
		{
			remaining_row_ind.insert(l);
		}
		for (int l = 0; l < N2; l++)
		{
			remaining_col_ind.insert(l);
		}
		int l_local;
		for (l_local = 0; l_local < N1; l_local++)
		{
			row_index = l_local;
			row = K->getRow(row_indices[row_index], col_indices);
			this->maxAbsVector(row, remaining_col_ind, max, col_index);
			if (abs(row(int(col_index))) > pow(10, -36))
			{
				break;
			}
		}
		if (l_local == N1)
		{
			return;
		}
		v = row;
		row_bases.push_back(row_index);
		col_bases.push_back(int(col_index));
		col = K->getCol(col_indices[col_index], row_indices);
		u = col / row(int(col_index));
		Uvec.push_back(u);
		Vvec.push_back(v);

		remaining_col_ind.erase(col_index);
		remaining_row_ind.erase(row_index);
		computed_rank = 1;

		double normS = 0.0;
		double prev_normS = __DBL_MAX__;
		this->maxAbsVector(col, remaining_row_ind, max, row_index);

		while (u.norm() * v.norm() > tol_ACA * normS && computed_rank < min)
		{
			row_bases.push_back(int(row_index));
			row = K->getRow(row_indices[row_index], col_indices);
			row_temp = row;
			for (int l = 0; l <= computed_rank - 1; ++l)
			{
				for (auto s = 0; s < Vvec[l].size(); s++)
				{
					row(s) -= Uvec[l](row_index) * Vvec[l](s);
				}
			}
			this->maxAbsVector(row, remaining_col_ind, max, col_index);
			col_bases.push_back(int(col_index));
			v = row;
			col = K->getCol(col_indices[col_index], row_indices);
			col_temp = col;
			for (int l = 0; l <= computed_rank - 1; ++l)
			{
				for (auto s = 0; s < Uvec[l].size(); s++)
				{
					col(s) -= Vvec[l](col_index) * Uvec[l](s);
				}
			}
			u = col / row(int(col_index));
			if (u.norm() < 1e-36 || v.norm() < 1e-36)
			{
				row_bases.pop_back();
				col_bases.pop_back();
				break;
			}
			Uvec.push_back(u);
			Vvec.push_back(v);

			++computed_rank;
			remaining_col_ind.erase(col_index);
			remaining_row_ind.erase(row_index);
			if (computed_rank != 2)
				prev_normS = normS;
			normS = pow(normS, 2) + pow(u.norm() * v.norm(), 2);
			for (int l = 0; l <= computed_rank - 1; ++l)
			{
				normS += 2 * abs(Uvec[l].dot(u) * Vvec[l].dot(v));
			}
			normS = sqrt(normS);
			this->maxAbsVector(col, remaining_row_ind, max, row_index);
		}

		L = Mat::Zero(computed_rank, computed_rank);
		R = Mat::Zero(computed_rank, computed_rank);
		if (computed_rank > 0)
		{
			for (auto i = 0; i < computed_rank; i++)
			{
				L(i, i) = 1.0;
				if (i >= 1)
				{
					for (auto j = 0; j <= i - 1; j++)
					{
						L(i, j) = Uvec[j](row_bases[i]);
					}
				}
			}
			for (auto i = 0; i < computed_rank; i++)
			{
				R(i, i) = Vvec[i](col_bases[i]);
				if (i >= 1)
				{
					for (auto j = 0; j <= i - 1; j++)
					{
						R(j, i) = Vvec[j](col_bases[i]);
					}
				}
			}
		}
	}

	void ACA_only_nodes(std::vector<int> &row_bases, std::vector<int> &col_bases, int &computed_rank, Mat &Ac, Mat &Ar, Mat &L, Mat &R)
	{
		int col_index;
		int row_index;
		int N1 = row_indices.size();
		int N2 = col_indices.size();
		Vec row(N2), col(N1), v(N2), u(N1);
		Vec row_temp, col_temp;
		std::vector<Vec> AcVec;
		std::vector<Vec> ArVec;

		std::vector<Vec> Uvec;
		std::vector<Vec> Vvec;
		computed_rank = 0;
		dtype max;
		int min = N1;
		if (N1 > N2)
		{
			min = N2;
		}

		std::set<int> remaining_row_ind;
		std::set<int> remaining_col_ind;

		for (int l = 0; l < N1; l++)
		{
			remaining_row_ind.insert(l);
		}
		for (int l = 0; l < N2; l++)
		{
			remaining_col_ind.insert(l);
		}
		int l_local;
		for (l_local = 0; l_local < N1; l_local++)
		{
			row_index = l_local;
			row = K->getRow(row_indices[row_index], col_indices);
			this->maxAbsVector(row, remaining_col_ind, max, col_index);
			if (abs(row(int(col_index))) > pow(10, -36))
			{
				break;
			}
		}
		if (l_local == N1)
		{
			Ac = Mat(N1, computed_rank);
			Ar = Mat(computed_rank, N2);
			return;
		}
		v = row;
		row_bases.push_back(row_index);
		col_bases.push_back(int(col_index));
		col = K->getCol(col_indices[col_index], row_indices);
		u = col / row(int(col_index));
		Uvec.push_back(u);
		Vvec.push_back(v);
		AcVec.push_back(col);
		ArVec.push_back(row);

		remaining_col_ind.erase(col_index);
		remaining_row_ind.erase(row_index);
		computed_rank = 1;

		double normS = 0.0;
		double prev_normS = __DBL_MAX__;
		this->maxAbsVector(col, remaining_row_ind, max, row_index);

		while (u.norm() * v.norm() > tol_ACA * normS && computed_rank < min)
		{
			row_bases.push_back(int(row_index));
			row = K->getRow(row_indices[row_index], col_indices);
			row_temp = row;
			for (int l = 0; l <= computed_rank - 1; ++l)
			{
				for (auto s = 0; s < Vvec[l].size(); s++)
				{
					row(s) -= Uvec[l](row_index) * Vvec[l](s);
				}
			}
			this->maxAbsVector(row, remaining_col_ind, max, col_index);
			col_bases.push_back(int(col_index));
			v = row;
			col = K->getCol(col_indices[col_index], row_indices);
			col_temp = col;
			for (int l = 0; l <= computed_rank - 1; ++l)
			{
				for (auto s = 0; s < Uvec[l].size(); s++)
				{
					col(s) -= Vvec[l](col_index) * Uvec[l](s);
				}
			}
			u = col / row(int(col_index));
			if (u.norm() < 1e-36 || v.norm() < 1e-36)
			{
				row_bases.pop_back();
				col_bases.pop_back();
				break;
			}
			Uvec.push_back(u);
			Vvec.push_back(v);
			AcVec.push_back(col_temp);
			ArVec.push_back(row_temp);

			++computed_rank;
			remaining_col_ind.erase(col_index);
			remaining_row_ind.erase(row_index);
			if (computed_rank != 2)
				prev_normS = normS;
			normS = pow(normS, 2) + pow(u.norm() * v.norm(), 2);
			for (int l = 0; l <= computed_rank - 1; ++l)
			{
				normS += 2 * abs(Uvec[l].dot(u) * Vvec[l].dot(v));
			}
			normS = sqrt(normS);
			this->maxAbsVector(col, remaining_row_ind, max, row_index);
		}
		Ac = Mat(N1, computed_rank);
		Ar = Mat(computed_rank, N2);
		for (auto i = 0; i < computed_rank; i++)
		{
			Ac.col(i) = AcVec[i];
			Ar.row(i) = ArVec[i];
		}
		L = Mat::Zero(computed_rank, computed_rank);
		R = Mat::Zero(computed_rank, computed_rank);
		if (computed_rank > 0)
		{
			for (auto i = 0; i < computed_rank; i++)
			{
				L(i, i) = 1.0;
				if (i >= 1)
				{
					for (auto j = 0; j <= i - 1; j++)
					{
						L(i, j) = Uvec[j](row_bases[i]);
					}
				}
			}
			for (auto i = 0; i < computed_rank; i++)
			{
				R(i, i) = Vvec[i](col_bases[i]);
				if (i >= 1)
				{
					for (auto j = 0; j <= i - 1; j++)
					{
						R(j, i) = Vvec[j](col_bases[i]);
					}
				}
			}
		}
	}
};
#endif
