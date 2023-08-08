//
// Remodified from : https://github.com/vaishna77/HODLR3D/ to keep comparable with FMM or H matrix
// HODLR3D Tree
//
#ifndef _FMM3DTreeRAMeff2_HPP__
#define _FMM3DTreeRAMeff2_HPP__
#include "kernel.hpp"
#include "aca.hpp"
#include "HODLR3Dnode.hpp"

template <typename kerneltype>
class HODLR3DTree
{
public:
	kerneltype *K;
	int nLevels;			//	Number of levels in the tree.
	int N;					//	Number of particles.
	int cubeRootN;			//	Number of particles.
	double L;				//	Semi-length of the simulation box.
	double smallestBoxSize; //	This is L/2.0^(nLevels).

	std::vector<int> nBoxesPerLevel;		   //	Number of boxes at each level in the tree.
	std::vector<double> boxRadius;			   //	Box radius at each level in the tree assuming the box at the root is [-1,1]^2
	std::vector<std::vector<HODLR3DNode*> > tree; //	The tree storing all the information.

	double ACA_epsilon;
	int nParticlesInLeafAlong1D;
	int nParticlesInLeaf;
	std::vector<double> Nodes1D;
	std::vector<pts3D> Nodes;
	std::vector<pts3D> gridPoints; // all particles in domain
	int TOL_POW;
	double elapsed_mvp = 0.0;
	double elapsed_assem = 0.0;
	double memory = 0.0;
	Vec collective_potential, collective_charge, computed_potential;

	// public:
	HODLR3DTree(kerneltype *K, int cubeRootN, int nLevels, int nParticlesInLeafAlong1D, double L, int TOL_POW)
	{
		this->K = K;
		this->cubeRootN = cubeRootN;
		this->nLevels = nLevels;
		this->L = L;
		this->nParticlesInLeafAlong1D = nParticlesInLeafAlong1D;
		this->nParticlesInLeaf = nParticlesInLeafAlong1D * nParticlesInLeafAlong1D * nParticlesInLeafAlong1D;
		this->TOL_POW = TOL_POW;
		nBoxesPerLevel.push_back(1);
		boxRadius.push_back(L);
		for (int k = 1; k <= nLevels; ++k)
		{
			nBoxesPerLevel.push_back(8 * nBoxesPerLevel[k - 1]);
			boxRadius.push_back(0.5 * boxRadius[k - 1]);
		}
		this->smallestBoxSize = boxRadius[nLevels];
		K->a = smallestBoxSize;
		this->N = cubeRootN * cubeRootN * cubeRootN;
	}

	void set_Standard_Cheb_Nodes()
	{
		for (int k = 0; k < cubeRootN; ++k)
		{
			Nodes1D.push_back(-cos((k + 0.5) / cubeRootN * PI));
		}
		pts3D temp1;
		for (int j = 0; j < cubeRootN; ++j)
		{
			for (int k = 0; k < cubeRootN; ++k)
			{
				for (int i = 0; i < cubeRootN; ++i)
				{
					temp1.x = Nodes1D[k];
					temp1.y = Nodes1D[j];
					temp1.z = Nodes1D[i];
					K->particles.push_back(temp1);
				}
			}
		}
	}

	void set_Uniform_Nodes()
	{
		for (int k = 0; k < cubeRootN; ++k)
		{
			Nodes1D.push_back(-L + 2.0 * L * (k + 1.0) / (cubeRootN + 1.0));
		}
		pts3D temp1;
		for (int j = 0; j < cubeRootN; ++j)
		{
			for (int k = 0; k < cubeRootN; ++k)
			{
				for (int i = 0; i < cubeRootN; ++i)
				{
					temp1.x = Nodes1D[k];
					temp1.y = Nodes1D[j];
					temp1.z = Nodes1D[i];
					K->particles.push_back(temp1);
				}
			}
		}
	}

	void shift_Nodes(double radius, pts3D center, std::vector<pts3D> &particle_loc)
	{
		for (size_t i = 0; i < Nodes.size(); i++)
		{
			pts3D temp;
			temp.x = Nodes[i].x * radius + center.x;
			temp.y = Nodes[i].y * radius + center.y;
			temp.z = Nodes[i].z * radius + center.z;
			particle_loc.push_back(temp);
		}
	}

	void createTree()
	{
		//	First create root and add to tree
		HODLR3DNode* root = new HODLR3DNode;
		root->boxNumber = 0;
		root->parentNumber = -1;
		#pragma omp parallel for
		for (int l = 0; l < 8; ++l)
		{
			root->childrenNumbers[l] = l;
		}
		#pragma omp parallel for
		for (int l = 0; l < 26; ++l)
		{
			root->neighborNumbers[l] = -1;
		}
		std::vector<HODLR3DNode*> rootLevel;
		rootLevel.push_back(root);
		tree.push_back(rootLevel);

		for (int j = 1; j <= nLevels; ++j)
		{
			std::vector<HODLR3DNode*> level;
			for (int k = 0; k < nBoxesPerLevel[j]; ++k)
			{
				HODLR3DNode* box = new HODLR3DNode;
				box->boxNumber = k;
				box->parentNumber = k / 8;
				for (int l = 0; l < 8; ++l)
				{
					box->childrenNumbers[l] = 8 * k + l;
				}
				level.push_back(box);
			}
			tree.push_back(level);
		}
	}

	//	Assigns the interactions for child0 of a box
	void assign_Child0_Interaction(int j, int k)
	{
		int nL = j + 1;
		int nC = 8 * k;
		int nN;

		//	Assign siblings
		{
			tree[nL][nC]->neighborNumbers[13] = 8 * k + 1;
			tree[nL][nC]->neighborNumbers[16] = 8 * k + 2;
			tree[nL][nC]->neighborNumbers[15] = 8 * k + 3;
			tree[nL][nC]->neighborNumbers[21] = 8 * k + 4;
			tree[nL][nC]->neighborNumbers[22] = 8 * k + 5;
			tree[nL][nC]->interactionList.push_back(8 * k + 6);
			tree[nL][nC]->neighborNumbers[24] = 8 * k + 7;
		}

		//	Assign children of parent's first neighbor
		{
			nN = tree[j][k]->neighborNumbers[1];
			if (nN != -1)
			{
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[6]);
				tree[nL][nC]->neighborNumbers[1] = tree[j][nN]->childrenNumbers[7];
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[0]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[1]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[2]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[3]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[4]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[5]);
			}
		}

		//	Assign children of parent's third neighbor
		{
			nN = tree[j][k]->neighborNumbers[3];
			if (nN != -1)
			{
				tree[nL][nC]->neighborNumbers[3] = tree[j][nN]->childrenNumbers[5];
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[6]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[0]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[1]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[2]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[3]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[4]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[7]);
			}
		}

		//	Assign children of parent's fourth neighbor
		{
			nN = tree[j][k]->neighborNumbers[4];
			if (nN != -1)
			{
				tree[nL][nC]->neighborNumbers[4] = tree[j][nN]->childrenNumbers[4];
				tree[nL][nC]->neighborNumbers[5] = tree[j][nN]->childrenNumbers[5];
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[6]);
				tree[nL][nC]->neighborNumbers[7] = tree[j][nN]->childrenNumbers[7];
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[0]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[1]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[2]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[3]);
			}
		}

		//	Assign children of parent's fifth neighbor
		{
			nN = tree[j][k]->neighborNumbers[5];
			if (nN != -1)
			{
				for (size_t i = 0; i < 8; i++)
				{
					tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[i]);
				}
			}
		}

		//	Assign children of parent's seventh neighbor
		{
			nN = tree[j][k]->neighborNumbers[7];
			if (nN != -1)
			{
				for (size_t i = 0; i < 8; i++)
				{
					tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[i]);
				}
			}
		}

		//	Assign children of parent's ninth neighbor
		{
			nN = tree[j][k]->neighborNumbers[9];
			if (nN != -1)
			{
				tree[nL][nC]->neighborNumbers[9] = tree[j][nN]->childrenNumbers[2];
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[6]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[0]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[1]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[3]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[4]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[5]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[7]);
			}
		}

		//	Assign children of parent's tenth neighbor
		{
			nN = tree[j][k]->neighborNumbers[10];
			if (nN != -1)
			{
				tree[nL][nC]->neighborNumbers[11] = tree[j][nN]->childrenNumbers[2];
				tree[nL][nC]->neighborNumbers[10] = tree[j][nN]->childrenNumbers[3];
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[6]);
				tree[nL][nC]->neighborNumbers[18] = tree[j][nN]->childrenNumbers[7];
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[0]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[1]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[4]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[5]);
			}
		}

		//	Assign children of parent's eleventh neighbor
		{
			nN = tree[j][k]->neighborNumbers[11];
			if (nN != -1)
			{
				for (size_t i = 0; i < 8; i++)
				{
					tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[i]);
				}
			}
		}

		//	Assign children of parent's twelveth neighbor
		{
			nN = tree[j][k]->neighborNumbers[12];
			if (nN != -1)
			{
				tree[nL][nC]->neighborNumbers[12] = tree[j][nN]->childrenNumbers[1];
				tree[nL][nC]->neighborNumbers[14] = tree[j][nN]->childrenNumbers[2];
				tree[nL][nC]->neighborNumbers[20] = tree[j][nN]->childrenNumbers[5];
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[6]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[0]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[3]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[4]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[7]);
			}
		}

		for (size_t n = 13; n <= 25; n++)
		{
			//	Assign children of parent's nth neighbor
			{
				if (n == 17 || n == 19 || n == 23 || n == 25)
				{
					continue;
				}
				else
				{
					nN = tree[j][k]->neighborNumbers[n];
					if (nN != -1)
					{
						for (size_t i = 0; i < 8; i++)
						{
							tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[i]);
						}
					}
				}
			}
		}
	}

	//	Assigns the interactions for child1 of a box
	void assign_Child1_Interaction(int j, int k)
	{
		int nL = j + 1;
		int nC = 8 * k + 1;
		int nN;

		//	Assign siblings
		{
			tree[nL][nC]->neighborNumbers[12] = 8 * k + 0;
			tree[nL][nC]->neighborNumbers[15] = 8 * k + 2;
			tree[nL][nC]->neighborNumbers[14] = 8 * k + 3;
			tree[nL][nC]->neighborNumbers[20] = 8 * k + 4;
			tree[nL][nC]->neighborNumbers[21] = 8 * k + 5;
			tree[nL][nC]->neighborNumbers[24] = 8 * k + 6;
			tree[nL][nC]->interactionList.push_back(8 * k + 7);
		}

		//	Assign children of parent's first neighbor
		{
			nN = tree[j][k]->neighborNumbers[1];
			if (nN != -1)
			{
				tree[nL][nC]->neighborNumbers[1] = tree[j][nN]->childrenNumbers[6];
				// tree[nL][nC]->neighborNumbers[0]	=	tree[j][nN]->childrenNumbers[7];
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[7]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[0]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[1]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[2]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[3]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[4]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[5]);
			}
		}

		//	Assign children of parent's third neighbor
		{
			nN = tree[j][k]->neighborNumbers[3];
			if (nN != -1)
			{
				for (size_t i = 0; i < 8; i++)
				{
					tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[i]);
				}
			}
		}

		//	Assign children of parent's fourth neighbor
		{
			nN = tree[j][k]->neighborNumbers[4];
			if (nN != -1)
			{
				tree[nL][nC]->neighborNumbers[3] = tree[j][nN]->childrenNumbers[4];
				tree[nL][nC]->neighborNumbers[4] = tree[j][nN]->childrenNumbers[5];
				tree[nL][nC]->neighborNumbers[7] = tree[j][nN]->childrenNumbers[6];
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[7]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[0]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[1]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[2]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[3]);
			}
		}

		//	Assign children of parent's fifth neighbor
		{
			nN = tree[j][k]->neighborNumbers[5];
			if (nN != -1)
			{
				tree[nL][nC]->neighborNumbers[5] = tree[j][nN]->childrenNumbers[4];
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[7]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[0]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[1]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[2]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[3]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[5]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[6]);
			}
		}

		//	Assign children of parent's seventh neighbor
		{
			nN = tree[j][k]->neighborNumbers[7];
			if (nN != -1)
			{
				for (size_t i = 0; i < 8; i++)
				{
					tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[i]);
				}
			}
		}

		//	Assign children of parent's ninth neighbor
		{
			nN = tree[j][k]->neighborNumbers[9];
			if (nN != -1)
			{
				for (size_t i = 0; i < 8; i++)
				{
					tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[i]);
				}
			}
		}

		//	Assign children of parent's tenth neighbor
		{
			nN = tree[j][k]->neighborNumbers[10];
			if (nN != -1)
			{
				tree[nL][nC]->neighborNumbers[10] = tree[j][nN]->childrenNumbers[2];
				tree[nL][nC]->neighborNumbers[9] = tree[j][nN]->childrenNumbers[3];
				tree[nL][nC]->neighborNumbers[18] = tree[j][nN]->childrenNumbers[6];
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[7]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[0]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[1]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[4]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[5]);
			}
		}

		//	Assign children of parent's eleventh neighbor
		{
			nN = tree[j][k]->neighborNumbers[11];
			if (nN != -1)
			{
				tree[nL][nC]->neighborNumbers[11] = tree[j][nN]->childrenNumbers[3];
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[7]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[0]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[1]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[2]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[4]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[5]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[6]);
			}
		}

		//	Assign children of parent's twelveth neighbor
		{
			nN = tree[j][k]->neighborNumbers[12];
			if (nN != -1)
			{
				for (size_t i = 0; i < 8; i++)
				{
					tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[i]);
				}
			}
		}

		//	Assign children of parent's thirteenth neighbor
		{
			nN = tree[j][k]->neighborNumbers[13];
			if (nN != -1)
			{
				tree[nL][nC]->neighborNumbers[13] = tree[j][nN]->childrenNumbers[0];
				tree[nL][nC]->neighborNumbers[16] = tree[j][nN]->childrenNumbers[3];
				tree[nL][nC]->neighborNumbers[22] = tree[j][nN]->childrenNumbers[4];
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[7]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[1]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[2]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[5]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[6]);
			}
		}

		for (size_t n = 14; n <= 25; n++)
		{
			//	Assign children of parent's nth neighbor
			{
				if (n == 17 || n == 19 || n == 23 || n == 25)
				{
					continue;
				}
				else
				{
					nN = tree[j][k]->neighborNumbers[n];
					if (nN != -1)
					{
						for (size_t i = 0; i < 8; i++)
						{
							tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[i]);
						}
					}
				}
			}
		}
	}

	//	Assigns the interactions for child2 of a box
	void assign_Child2_Interaction(int j, int k)
	{
		int nL = j + 1;
		int nC = 8 * k + 2;
		int nN;

		//	Assign siblings
		{
			tree[nL][nC]->neighborNumbers[9] = 8 * k + 0;
			tree[nL][nC]->neighborNumbers[10] = 8 * k + 1;
			tree[nL][nC]->neighborNumbers[12] = 8 * k + 3;
			tree[nL][nC]->interactionList.push_back(8 * k + 4);
			tree[nL][nC]->neighborNumbers[18] = 8 * k + 5;
			tree[nL][nC]->neighborNumbers[21] = 8 * k + 6;
			tree[nL][nC]->neighborNumbers[20] = 8 * k + 7;
		}

		//	Assign children of parent's first neighbor
		{
			nN = tree[j][k]->neighborNumbers[1];
			if (nN != -1)
			{
				for (size_t i = 0; i < 8; i++)
				{
					tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[i]);
				}
			}
		}

		//	Assign children of parent's third neighbor
		{
			nN = tree[j][k]->neighborNumbers[3];
			if (nN != -1)
			{
				for (size_t i = 0; i < 8; i++)
				{
					tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[i]);
				}
			}
		}

		//	Assign children of parent's fourth neighbor
		{
			nN = tree[j][k]->neighborNumbers[4];
			if (nN != -1)
			{
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[4]);
				tree[nL][nC]->neighborNumbers[1] = tree[j][nN]->childrenNumbers[5];
				tree[nL][nC]->neighborNumbers[4] = tree[j][nN]->childrenNumbers[6];
				tree[nL][nC]->neighborNumbers[3] = tree[j][nN]->childrenNumbers[7];
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[0]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[1]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[2]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[3]);
			}
		}

		//	Assign children of parent's fifth neighbor
		{
			nN = tree[j][k]->neighborNumbers[5];
			if (nN != -1)
			{
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[4]);
				tree[nL][nC]->neighborNumbers[5] = tree[j][nN]->childrenNumbers[7];
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[0]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[1]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[2]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[3]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[5]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[6]);
			}
		}

		//	Assign children of parent's seventh neighbor
		{
			nN = tree[j][k]->neighborNumbers[7];
			if (nN != -1)
			{
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[4]);
				tree[nL][nC]->neighborNumbers[7] = tree[j][nN]->childrenNumbers[5];
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[0]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[1]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[2]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[3]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[6]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[7]);
			}
		}

		//	Assign children of parent's ninth neighbor
		{
			nN = tree[j][k]->neighborNumbers[9];
			if (nN != -1)
			{
				for (size_t i = 0; i < 8; i++)
				{
					tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[i]);
				}
			}
		}

		//	Assign children of parent's tenth neighbor
		{
			nN = tree[j][k]->neighborNumbers[10];
			if (nN != -1)
			{
				for (size_t i = 0; i < 8; i++)
				{
					tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[i]);
				}
			}
		}

		//	Assign children of parent's eleventh neighbor
		{
			nN = tree[j][k]->neighborNumbers[11];
			if (nN != -1)
			{
				for (size_t i = 0; i < 8; i++)
				{
					tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[i]);
				}
			}
		}

		//	Assign children of parent's twelveth neighbor
		{
			nN = tree[j][k]->neighborNumbers[12];
			if (nN != -1)
			{
				for (size_t i = 0; i < 8; i++)
				{
					tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[i]);
				}
			}
		}

		//	Assign children of parent's thirteenth neighbor
		{
			nN = tree[j][k]->neighborNumbers[13];
			if (nN != -1)
			{
				tree[nL][nC]->neighborNumbers[11] = tree[j][nN]->childrenNumbers[0];
				tree[nL][nC]->neighborNumbers[13] = tree[j][nN]->childrenNumbers[3];
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[4]);
				tree[nL][nC]->neighborNumbers[22] = tree[j][nN]->childrenNumbers[7];
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[1]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[2]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[5]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[6]);
			}
		}

		//	Assign children of parent's fourteenth neighbor
		{
			nN = tree[j][k]->neighborNumbers[14];
			if (nN != -1)
			{
				for (size_t i = 0; i < 8; i++)
				{
					tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[i]);
				}
			}
		}

		//	Assign children of parent's fifteenth neighbor
		{
			nN = tree[j][k]->neighborNumbers[15];
			if (nN != -1)
			{
				tree[nL][nC]->neighborNumbers[14] = tree[j][nN]->childrenNumbers[0];
				tree[nL][nC]->neighborNumbers[15] = tree[j][nN]->childrenNumbers[1];
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[4]);
				tree[nL][nC]->neighborNumbers[24] = tree[j][nN]->childrenNumbers[5];
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[2]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[3]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[6]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[7]);
			}
		}

		//	Assign children of parent's sixteenth neighbor
		{
			nN = tree[j][k]->neighborNumbers[16];
			if (nN != -1)
			{
				tree[nL][nC]->neighborNumbers[16] = tree[j][nN]->childrenNumbers[0];
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[4]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[1]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[2]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[3]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[5]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[6]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[7]);
			}
		}

		for (size_t n = 17; n <= 25; n++)
		{
			//	Assign children of parent's nth neighbor
			{
				if (n == 17 || n == 19 || n == 23 || n == 25)
				{
					continue;
				}
				else
				{
					nN = tree[j][k]->neighborNumbers[n];
					if (nN != -1)
					{
						for (size_t i = 0; i < 8; i++)
						{
							tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[i]);
						}
					}
				}
			}
		}
	}

	//	Assigns the interactions for child3 of a box
	void assign_Child3_Interaction(int j, int k)
	{
		int nL = j + 1;
		int nC = 8 * k + 3;
		int nN;

		//	Assign siblings
		{
			tree[nL][nC]->neighborNumbers[10] = 8 * k + 0;
			tree[nL][nC]->neighborNumbers[11] = 8 * k + 1;
			tree[nL][nC]->neighborNumbers[13] = 8 * k + 2;
			tree[nL][nC]->neighborNumbers[18] = 8 * k + 4;
			tree[nL][nC]->interactionList.push_back(8 * k + 5);
			tree[nL][nC]->neighborNumbers[22] = 8 * k + 6;
			tree[nL][nC]->neighborNumbers[21] = 8 * k + 7;
		}

		for (size_t n = 1; n <= 1; n++)
		{
			//	Assign children of parent's nth neighbor
			{
				nN = tree[j][k]->neighborNumbers[n];
				if (nN != -1)
				{
					for (size_t i = 0; i < 8; i++)
					{
						tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[i]);
					}
				}
			}
		}

		//	Assign children of parent's third neighbor
		{
			nN = tree[j][k]->neighborNumbers[3];
			if (nN != -1)
			{
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[5]);
				tree[nL][nC]->neighborNumbers[3] = tree[j][nN]->childrenNumbers[6];
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[0]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[1]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[2]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[3]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[4]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[7]);
			}
		}

		//	Assign children of parent's fourth neighbor
		{
			nN = tree[j][k]->neighborNumbers[4];
			if (nN != -1)
			{
				tree[nL][nC]->neighborNumbers[1] = tree[j][nN]->childrenNumbers[4];
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[5]);
				tree[nL][nC]->neighborNumbers[5] = tree[j][nN]->childrenNumbers[6];
				tree[nL][nC]->neighborNumbers[4] = tree[j][nN]->childrenNumbers[7];
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[0]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[1]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[2]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[3]);
			}
		}

		//	Assign children of parent's fifth neighbor
		{
			nN = tree[j][k]->neighborNumbers[5];
			if (nN != -1)
			{
				for (size_t i = 0; i < 8; i++)
				{
					tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[i]);
				}
			}
		}

		//	Assign children of parent's seventh neighbor
		{
			nN = tree[j][k]->neighborNumbers[7];
			if (nN != -1)
			{
				tree[nL][nC]->neighborNumbers[7] = tree[j][nN]->childrenNumbers[4];
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[5]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[0]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[1]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[2]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[3]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[6]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[7]);
			}
		}

		for (size_t n = 9; n <= 11; n++)
		{
			//	Assign children of parent's nth neighbor
			{
				nN = tree[j][k]->neighborNumbers[n];
				if (nN != -1)
				{
					for (size_t i = 0; i < 8; i++)
					{
						tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[i]);
					}
				}
			}
		}

		//	Assign children of parent's twelveth neighbor
		{
			nN = tree[j][k]->neighborNumbers[12];
			if (nN != -1)
			{
				tree[nL][nC]->neighborNumbers[9] = tree[j][nN]->childrenNumbers[1];
				tree[nL][nC]->neighborNumbers[12] = tree[j][nN]->childrenNumbers[2];
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[5]);
				tree[nL][nC]->neighborNumbers[20] = tree[j][nN]->childrenNumbers[6];
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[0]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[3]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[4]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[7]);
			}
		}

		//	Assign children of parent's thirteenth neighbor
		{
			nN = tree[j][k]->neighborNumbers[13];
			if (nN != -1)
			{
				for (size_t i = 0; i < 8; i++)
				{
					tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[i]);
				}
			}
		}

		//	Assign children of parent's fourteenth neighbor
		{
			nN = tree[j][k]->neighborNumbers[14];
			if (nN != -1)
			{
				tree[nL][nC]->neighborNumbers[14] = tree[j][nN]->childrenNumbers[1];
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[5]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[0]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[2]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[3]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[4]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[6]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[7]);
			}
		}

		//	Assign children of parent's fifteenth neighbor
		{
			nN = tree[j][k]->neighborNumbers[15];
			if (nN != -1)
			{
				tree[nL][nC]->neighborNumbers[15] = tree[j][nN]->childrenNumbers[0];
				tree[nL][nC]->neighborNumbers[16] = tree[j][nN]->childrenNumbers[1];
				tree[nL][nC]->neighborNumbers[24] = tree[j][nN]->childrenNumbers[4];
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[5]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[2]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[3]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[6]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[7]);
			}
		}

		for (size_t n = 16; n <= 25; n++)
		{
			//	Assign children of parent's nth neighbor
			{
				if (n == 17 || n == 19 || n == 23 || n == 25)
				{
					continue;
				}
				else
				{
					nN = tree[j][k]->neighborNumbers[n];
					if (nN != -1)
					{
						for (size_t i = 0; i < 8; i++)
						{
							tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[i]);
						}
					}
				}
			}
		}
	}

	//	Assigns the interactions for child4 of a box
	void assign_Child4_Interaction(int j, int k)
	{
		int nL = j + 1;
		int nC = 8 * k + 4;
		int nN;

		//	Assign siblings
		{
			tree[nL][nC]->neighborNumbers[4] = 8 * k + 0;
			tree[nL][nC]->neighborNumbers[5] = 8 * k + 1;
			tree[nL][nC]->interactionList.push_back(8 * k + 2);
			tree[nL][nC]->neighborNumbers[7] = 8 * k + 3;
			tree[nL][nC]->neighborNumbers[13] = 8 * k + 5;
			tree[nL][nC]->neighborNumbers[16] = 8 * k + 6;
			tree[nL][nC]->neighborNumbers[15] = 8 * k + 7;
		}

		for (size_t n = 1; n <= 8; n++)
		{
			//	Assign children of parent's nth neighbor
			{
				if (n == 2 || n == 6 || n == 8)
				{
					continue;
				}
				else
				{
					nN = tree[j][k]->neighborNumbers[n];
					if (nN != -1)
					{
						for (size_t i = 0; i < 8; i++)
						{
							tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[i]);
						}
					}
				}
			}
		}

		//	Assign children of parent's nineth neighbor
		{
			nN = tree[j][k]->neighborNumbers[9];
			if (nN != -1)
			{
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[2]);
				tree[nL][nC]->neighborNumbers[9] = tree[j][nN]->childrenNumbers[6];
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[0]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[1]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[3]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[4]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[5]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[7]);
			}
		}

		//	Assign children of parent's tenth neighbor
		{
			nN = tree[j][k]->neighborNumbers[10];
			if (nN != -1)
			{
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[2]);
				tree[nL][nC]->neighborNumbers[1] = tree[j][nN]->childrenNumbers[3];
				tree[nL][nC]->neighborNumbers[11] = tree[j][nN]->childrenNumbers[6];
				tree[nL][nC]->neighborNumbers[10] = tree[j][nN]->childrenNumbers[7];
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[0]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[1]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[4]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[5]);
			}
		}

		//	Assign children of parent's eleventh neighbor
		{
			nN = tree[j][k]->neighborNumbers[11];
			if (nN != -1)
			{
				for (size_t i = 0; i < 8; i++)
				{
					tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[i]);
				}
			}
		}

		//	Assign children of parent's twelveth neighbor
		{
			nN = tree[j][k]->neighborNumbers[12];
			if (nN != -1)
			{
				tree[nL][nC]->neighborNumbers[3] = tree[j][nN]->childrenNumbers[1];
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[2]);
				tree[nL][nC]->neighborNumbers[12] = tree[j][nN]->childrenNumbers[5];
				tree[nL][nC]->neighborNumbers[14] = tree[j][nN]->childrenNumbers[6];
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[0]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[3]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[4]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[7]);
			}
		}

		for (size_t n = 13; n <= 16; n++)
		{
			//	Assign children of parent's nth neighbor
			{
				nN = tree[j][k]->neighborNumbers[n];
				if (nN != -1)
				{
					for (size_t i = 0; i < 8; i++)
					{
						tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[i]);
					}
				}
			}
		}

		//	Assign children of parent's eighteenth neighbor
		{
			nN = tree[j][k]->neighborNumbers[18];
			if (nN != -1)
			{
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[2]);
				tree[nL][nC]->neighborNumbers[18] = tree[j][nN]->childrenNumbers[3];
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[0]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[1]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[4]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[5]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[6]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[7]);
			}
		}

		//	Assign children of parent's twentieth neighbor
		{
			nN = tree[j][k]->neighborNumbers[20];
			if (nN != -1)
			{
				tree[nL][nC]->neighborNumbers[20] = tree[j][nN]->childrenNumbers[1];
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[2]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[0]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[3]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[4]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[5]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[6]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[7]);
			}
		}

		//	Assign children of parent's twentyfirst neighbor
		{
			nN = tree[j][k]->neighborNumbers[21];
			if (nN != -1)
			{
				tree[nL][nC]->neighborNumbers[21] = tree[j][nN]->childrenNumbers[0];
				tree[nL][nC]->neighborNumbers[22] = tree[j][nN]->childrenNumbers[1];
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[2]);
				tree[nL][nC]->neighborNumbers[24] = tree[j][nN]->childrenNumbers[3];
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[4]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[5]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[6]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[7]);
			}
		}

		for (size_t n = 22; n <= 25; n++)
		{
			//	Assign children of parent's nth neighbor
			{
				if (n == 23 || n == 25)
				{
					continue;
				}
				else
				{
					nN = tree[j][k]->neighborNumbers[n];
					if (nN != -1)
					{
						for (size_t i = 0; i < 8; i++)
						{
							tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[i]);
						}
					}
				}
			}
		}
	}

	//	Assigns the interactions for child5 of a box
	void assign_Child5_Interaction(int j, int k)
	{
		int nL = j + 1;
		int nC = 8 * k + 5;
		int nN;

		//	Assign siblings
		{
			tree[nL][nC]->neighborNumbers[3] = 8 * k + 0;
			tree[nL][nC]->neighborNumbers[4] = 8 * k + 1;
			tree[nL][nC]->neighborNumbers[7] = 8 * k + 2;
			tree[nL][nC]->interactionList.push_back(8 * k + 3);
			tree[nL][nC]->neighborNumbers[12] = 8 * k + 4;
			tree[nL][nC]->neighborNumbers[15] = 8 * k + 6;
			tree[nL][nC]->neighborNumbers[14] = 8 * k + 7;
		}

		for (size_t n = 1; n <= 8; n++)
		{
			//	Assign children of parent's nth neighbor
			{
				if (n == 2 || n == 6 || n == 8)
				{
					continue;
				}
				else
				{
					nN = tree[j][k]->neighborNumbers[n];
					if (nN != -1)
					{
						for (size_t i = 0; i < 8; i++)
						{
							tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[i]);
						}
					}
				}
			}
		}

		//	Assign children of parent's nineth neighbor
		{
			nN = tree[j][k]->neighborNumbers[9];
			if (nN != -1)
			{
				for (size_t i = 0; i < 8; i++)
				{
					tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[i]);
				}
			}
		}

		//	Assign children of parent's tenth neighbor
		{
			nN = tree[j][k]->neighborNumbers[10];
			if (nN != -1)
			{
				tree[nL][nC]->neighborNumbers[1] = tree[j][nN]->childrenNumbers[2];
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[3]);
				tree[nL][nC]->neighborNumbers[10] = tree[j][nN]->childrenNumbers[6];
				tree[nL][nC]->neighborNumbers[9] = tree[j][nN]->childrenNumbers[7];
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[0]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[1]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[4]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[5]);
			}
		}

		//	Assign children of parent's eleventh neighbor
		{
			nN = tree[j][k]->neighborNumbers[11];
			if (nN != -1)
			{
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[3]);
				tree[nL][nC]->neighborNumbers[11] = tree[j][nN]->childrenNumbers[7];
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[0]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[1]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[2]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[4]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[5]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[6]);
			}
		}

		//	Assign children of parent's 12th neighbor
		{
			nN = tree[j][k]->neighborNumbers[12];
			if (nN != -1)
			{
				for (size_t i = 0; i < 8; i++)
				{
					tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[i]);
				}
			}
		}

		//	Assign children of parent's 13 neighbor
		{
			nN = tree[j][k]->neighborNumbers[13];
			if (nN != -1)
			{
				tree[nL][nC]->neighborNumbers[5] = tree[j][nN]->childrenNumbers[0];
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[3]);
				tree[nL][nC]->neighborNumbers[13] = tree[j][nN]->childrenNumbers[4];
				tree[nL][nC]->neighborNumbers[16] = tree[j][nN]->childrenNumbers[7];
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[1]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[2]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[5]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[6]);
			}
		}

		for (size_t n = 14; n <= 16; n++)
		{
			//	Assign children of parent's nth neighbor
			{
				nN = tree[j][k]->neighborNumbers[n];
				if (nN != -1)
				{
					for (size_t i = 0; i < 8; i++)
					{
						tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[i]);
					}
				}
			}
		}

		//	Assign children of parent's 18th neighbor
		{
			nN = tree[j][k]->neighborNumbers[18];
			if (nN != -1)
			{
				tree[nL][nC]->neighborNumbers[18] = tree[j][nN]->childrenNumbers[2];
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[3]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[0]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[1]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[4]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[5]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[6]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[7]);
			}
		}

		//	Assign children of parent's 20th neighbor
		{
			nN = tree[j][k]->neighborNumbers[20];
			if (nN != -1)
			{
				for (size_t i = 0; i < 8; i++)
				{
					tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[i]);
				}
			}
		}

		//	Assign children of parent's 21st neighbor
		{
			nN = tree[j][k]->neighborNumbers[21];
			if (nN != -1)
			{
				tree[nL][nC]->neighborNumbers[20] = tree[j][nN]->childrenNumbers[0];
				tree[nL][nC]->neighborNumbers[21] = tree[j][nN]->childrenNumbers[1];
				tree[nL][nC]->neighborNumbers[24] = tree[j][nN]->childrenNumbers[2];
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[3]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[4]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[5]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[6]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[7]);
			}
		}

		//	Assign children of parent's 22nd neighbor
		{
			nN = tree[j][k]->neighborNumbers[22];
			if (nN != -1)
			{
				tree[nL][nC]->neighborNumbers[22] = tree[j][nN]->childrenNumbers[0];
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[3]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[1]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[2]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[4]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[5]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[6]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[7]);
			}
		}

		for (size_t n = 24; n <= 24; n++)
		{
			//	Assign children of parent's nth neighbor
			{
				nN = tree[j][k]->neighborNumbers[n];
				if (nN != -1)
				{
					for (size_t i = 0; i < 8; i++)
					{
						tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[i]);
					}
				}
			}
		}
	}

	//	Assigns the interactions for child6 of a box
	void assign_Child6_Interaction(int j, int k)
	{
		int nL = j + 1;
		int nC = 8 * k + 6;
		int nN;

		//	Assign siblings
		{
			tree[nL][nC]->interactionList.push_back(8 * k + 0);
			tree[nL][nC]->neighborNumbers[1] = 8 * k + 1;
			tree[nL][nC]->neighborNumbers[4] = 8 * k + 2;
			tree[nL][nC]->neighborNumbers[3] = 8 * k + 3;
			tree[nL][nC]->neighborNumbers[9] = 8 * k + 4;
			tree[nL][nC]->neighborNumbers[10] = 8 * k + 5;
			tree[nL][nC]->neighborNumbers[12] = 8 * k + 7;
		}

		for (size_t n = 1; n <= 12; n++)
		{
			//	Assign children of parent's nth neighbor
			{
				if (n == 2 || n == 6 || n == 8)
				{
					continue;
				}
				else
				{
					nN = tree[j][k]->neighborNumbers[n];
					if (nN != -1)
					{
						for (size_t i = 0; i < 8; i++)
						{
							tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[i]);
						}
					}
				}
			}
		}

		//	Assign children of parent's 13th neighbor
		{
			nN = tree[j][k]->neighborNumbers[13];
			if (nN != -1)
			{
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[0]);
				tree[nL][nC]->neighborNumbers[5] = tree[j][nN]->childrenNumbers[3];
				tree[nL][nC]->neighborNumbers[11] = tree[j][nN]->childrenNumbers[4];
				tree[nL][nC]->neighborNumbers[13] = tree[j][nN]->childrenNumbers[7];
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[1]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[2]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[5]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[6]);
			}
		}

		//	Assign children of parent's 14th neighbor
		{
			nN = tree[j][k]->neighborNumbers[14];
			if (nN != -1)
			{
				for (size_t i = 0; i < 8; i++)
				{
					tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[i]);
				}
			}
		}

		//	Assign children of parent's 15th neighbor
		{
			nN = tree[j][k]->neighborNumbers[15];
			if (nN != -1)
			{
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[0]);
				tree[nL][nC]->neighborNumbers[7] = tree[j][nN]->childrenNumbers[1];
				tree[nL][nC]->neighborNumbers[14] = tree[j][nN]->childrenNumbers[4];
				tree[nL][nC]->neighborNumbers[15] = tree[j][nN]->childrenNumbers[5];
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[2]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[3]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[6]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[7]);
			}
		}

		//	Assign children of parent's 16th neighbor
		{
			nN = tree[j][k]->neighborNumbers[16];
			if (nN != -1)
			{
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[0]);
				tree[nL][nC]->neighborNumbers[16] = tree[j][nN]->childrenNumbers[4];
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[1]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[2]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[3]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[5]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[6]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[7]);
			}
		}

		for (size_t n = 18; n <= 20; n++)
		{
			//	Assign children of parent's nth neighbor
			{
				if (n == 19)
				{
					continue;
				}
				else
				{
					nN = tree[j][k]->neighborNumbers[n];
					if (nN != -1)
					{
						for (size_t i = 0; i < 8; i++)
						{
							tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[i]);
						}
					}
				}
			}
		}

		//	Assign children of parent's 21st neighbor
		{
			nN = tree[j][k]->neighborNumbers[21];
			if (nN != -1)
			{
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[0]);
				tree[nL][nC]->neighborNumbers[18] = tree[j][nN]->childrenNumbers[1];
				tree[nL][nC]->neighborNumbers[21] = tree[j][nN]->childrenNumbers[2];
				tree[nL][nC]->neighborNumbers[20] = tree[j][nN]->childrenNumbers[3];
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[4]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[5]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[6]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[7]);
			}
		}

		//	Assign children of parent's 22nd neighbor
		{
			nN = tree[j][k]->neighborNumbers[22];
			if (nN != -1)
			{
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[0]);
				tree[nL][nC]->neighborNumbers[22] = tree[j][nN]->childrenNumbers[3];
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[1]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[2]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[4]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[5]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[6]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[7]);
			}
		}

		//	Assign children of parent's 24th neighbor
		{
			nN = tree[j][k]->neighborNumbers[24];
			if (nN != -1)
			{
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[0]);
				tree[nL][nC]->neighborNumbers[24] = tree[j][nN]->childrenNumbers[1];
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[2]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[3]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[4]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[5]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[6]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[7]);
			}
		}
	}

	//	Assigns the interactions for child7 of a box
	void assign_Child7_Interaction(int j, int k)
	{
		int nL = j + 1;
		int nC = 8 * k + 7;
		int nN;

		//	Assign siblings
		{
			tree[nL][nC]->neighborNumbers[1] = 8 * k + 0;
			tree[nL][nC]->interactionList.push_back(8 * k + 1);
			tree[nL][nC]->neighborNumbers[5] = 8 * k + 2;
			tree[nL][nC]->neighborNumbers[4] = 8 * k + 3;
			tree[nL][nC]->neighborNumbers[10] = 8 * k + 4;
			tree[nL][nC]->neighborNumbers[11] = 8 * k + 5;
			tree[nL][nC]->neighborNumbers[13] = 8 * k + 6;
		}

		for (size_t n = 1; n <= 11; n++)
		{
			//	Assign children of parent's nth neighbor
			{
				if (n == 2 || n == 6 || n == 8)
				{
					continue;
				}
				else
				{
					nN = tree[j][k]->neighborNumbers[n];
					if (nN != -1)
					{
						for (size_t i = 0; i < 8; i++)
						{
							tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[i]);
						}
					}
				}
			}
		}

		//	Assign children of parent's 12th neighbor
		{
			nN = tree[j][k]->neighborNumbers[12];
			if (nN != -1)
			{
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[1]);
				tree[nL][nC]->neighborNumbers[3] = tree[j][nN]->childrenNumbers[2];
				tree[nL][nC]->neighborNumbers[9] = tree[j][nN]->childrenNumbers[5];
				tree[nL][nC]->neighborNumbers[12] = tree[j][nN]->childrenNumbers[6];
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[0]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[3]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[4]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[7]);
			}
		}

		//	Assign children of parent's 13th neighbor
		{
			nN = tree[j][k]->neighborNumbers[13];
			if (nN != -1)
			{
				for (size_t i = 0; i < 8; i++)
				{
					tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[i]);
				}
			}
		}

		//	Assign children of parent's 14th neighbor
		{
			nN = tree[j][k]->neighborNumbers[14];
			if (nN != -1)
			{
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[1]);
				tree[nL][nC]->neighborNumbers[14] = tree[j][nN]->childrenNumbers[5];
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[0]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[2]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[3]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[4]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[6]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[7]);
			}
		}

		//	Assign children of parent's 15th neighbor
		{
			nN = tree[j][k]->neighborNumbers[15];
			if (nN != -1)
			{
				tree[nL][nC]->neighborNumbers[7] = tree[j][nN]->childrenNumbers[0];
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[1]);
				tree[nL][nC]->neighborNumbers[15] = tree[j][nN]->childrenNumbers[4];
				tree[nL][nC]->neighborNumbers[16] = tree[j][nN]->childrenNumbers[5];
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[2]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[3]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[6]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[7]);
			}
		}

		for (size_t n = 16; n <= 18; n++)
		{
			//	Assign children of parent's nth neighbor
			{
				if (n == 17)
				{
					continue;
				}
				else
				{
					nN = tree[j][k]->neighborNumbers[n];
					if (nN != -1)
					{
						for (size_t i = 0; i < 8; i++)
						{
							tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[i]);
						}
					}
				}
			}
		}

		//	Assign children of parent's 20th neighbor
		{
			nN = tree[j][k]->neighborNumbers[20];
			if (nN != -1)
			{
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[1]);
				tree[nL][nC]->neighborNumbers[20] = tree[j][nN]->childrenNumbers[2];
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[0]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[3]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[4]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[5]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[6]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[7]);
			}
		}

		//	Assign children of parent's 21st neighbor
		{
			nN = tree[j][k]->neighborNumbers[21];
			if (nN != -1)
			{
				tree[nL][nC]->neighborNumbers[18] = tree[j][nN]->childrenNumbers[0];
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[1]);
				tree[nL][nC]->neighborNumbers[22] = tree[j][nN]->childrenNumbers[2];
				tree[nL][nC]->neighborNumbers[21] = tree[j][nN]->childrenNumbers[3];
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[4]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[5]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[6]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[7]);
			}
		}

		//	Assign children of parent's 22nd neighbor
		{
			nN = tree[j][k]->neighborNumbers[22];
			if (nN != -1)
			{
				for (size_t i = 0; i < 8; i++)
				{
					tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[i]);
				}
			}
		}

		//	Assign children of parent's 24th neighbor
		{
			nN = tree[j][k]->neighborNumbers[24];
			if (nN != -1)
			{
				tree[nL][nC]->neighborNumbers[24] = tree[j][nN]->childrenNumbers[0];
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[1]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[2]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[3]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[4]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[5]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[6]);
				tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[7]);
			}
		}
	}

	//	Assigns the interactions for the children of a box
	void assign_Box_Interactions(int j, int k)
	{
		assign_Child0_Interaction(j, k);
		assign_Child1_Interaction(j, k);
		assign_Child2_Interaction(j, k);
		assign_Child3_Interaction(j, k);
		assign_Child4_Interaction(j, k);
		assign_Child5_Interaction(j, k);
		assign_Child6_Interaction(j, k);
		assign_Child7_Interaction(j, k);
	}

	//	Assigns the interactions for the children all boxes at a given level
	void assign_Level_Interactions(int j)
	{
		#pragma omp parallel for
		for (int k = 0; k < nBoxesPerLevel[j]; ++k)
		{
			assign_Box_Interactions(j, k);
		}
	}

	//	Assigns the interactions for the children all boxes in the tree
	void assign_Tree_Interactions()
	{
		for (int j = 0; j < nLevels; ++j)
		{
			assign_Level_Interactions(j);
		}
	}

	void assign_Center_Location()
	{
		int J;
		tree[0][0]->center.x = 0.0;
		tree[0][0]->center.y = 0.0;
		tree[0][0]->center.z = 0.0;
		for (int j = 0; j < nLevels; ++j)
		{
			J = j + 1;
			double shift = 0.5 * boxRadius[j];
			#pragma omp parallel for
			for (int k = 0; k < nBoxesPerLevel[j]; ++k)
			{
				tree[J][8 * k]->center.x = tree[j][k]->center.x - shift;
				tree[J][8 * k + 1]->center.x = tree[j][k]->center.x + shift;
				tree[J][8 * k + 2]->center.x = tree[j][k]->center.x + shift;
				tree[J][8 * k + 3]->center.x = tree[j][k]->center.x - shift;
				tree[J][8 * k + 4]->center.x = tree[j][k]->center.x - shift;
				tree[J][8 * k + 5]->center.x = tree[j][k]->center.x + shift;
				tree[J][8 * k + 6]->center.x = tree[j][k]->center.x + shift;
				tree[J][8 * k + 7]->center.x = tree[j][k]->center.x - shift;

				tree[J][8 * k]->center.y = tree[j][k]->center.y - shift;
				tree[J][8 * k + 1]->center.y = tree[j][k]->center.y - shift;
				tree[J][8 * k + 2]->center.y = tree[j][k]->center.y + shift;
				tree[J][8 * k + 3]->center.y = tree[j][k]->center.y + shift;
				tree[J][8 * k + 4]->center.y = tree[j][k]->center.y - shift;
				tree[J][8 * k + 5]->center.y = tree[j][k]->center.y - shift;
				tree[J][8 * k + 6]->center.y = tree[j][k]->center.y + shift;
				tree[J][8 * k + 7]->center.y = tree[j][k]->center.y + shift;

				tree[J][8 * k]->center.z = tree[j][k]->center.z - shift;
				tree[J][8 * k + 1]->center.z = tree[j][k]->center.z - shift;
				tree[J][8 * k + 2]->center.z = tree[j][k]->center.z - shift;
				tree[J][8 * k + 3]->center.z = tree[j][k]->center.z - shift;
				tree[J][8 * k + 4]->center.z = tree[j][k]->center.z + shift;
				tree[J][8 * k + 5]->center.z = tree[j][k]->center.z + shift;
				tree[J][8 * k + 6]->center.z = tree[j][k]->center.z + shift;
				tree[J][8 * k + 7]->center.z = tree[j][k]->center.z + shift;
			}
		}
	}

	// Check the vertex sharing interactions
	void separate_vertex_interaction()
	{
		for (int j = 1; j <= nLevels; j++)
		{
			for (int k = 0; k < nBoxesPerLevel[j]; k++)
			{
				for (size_t i = 0; i < tree[j][k]->interactionList.size(); i++)
				{
					int ki = tree[j][k]->interactionList[i];

					double temp = dist3D(tree[j][k]->center, tree[j][ki]->center);

					if (temp - (2.0 * sqrt(3) * boxRadius[j]) <= 1.0e-16)
					{
						tree[j][k]->isVertex[ki] = true;
					}
					else
					{
						tree[j][k]->isVertex[ki] = false;
					}
				}
			}
		}
	}

	void assignChargeLocations()
	{
		for (int i = 0; i < N; i++)
		{
			tree[0][0]->chargeLocations.push_back(i);
		}
		for (int j = 0; j < nLevels; j++)
		{ // assign particles to its children
			for (int k = 0; k < nBoxesPerLevel[j]; k++)
			{
				int J = j + 1;
				int Kp = 8 * k;
				for (size_t i = 0; i < tree[j][k]->chargeLocations.size(); i++)
				{
					int index = tree[j][k]->chargeLocations[i];
					if (K->particles[index].z <= tree[j][k]->center.z)
					{ // children 0,1,2,3
						if (K->particles[index].x <= tree[j][k]->center.x)
						{ // children 0,3
							if (K->particles[index].y <= tree[j][k]->center.y)
							{ // child 0
								tree[J][Kp]->chargeLocations.push_back(index);
							}
							else
							{ // child 3
								tree[J][Kp + 3]->chargeLocations.push_back(index);
							}
						}
						else
						{ // children 1,2
							if (K->particles[index].y <= tree[j][k]->center.y)
							{ // child 1
								tree[J][Kp + 1]->chargeLocations.push_back(index);
							}
							else
							{ // child 2
								tree[J][Kp + 2]->chargeLocations.push_back(index);
							}
						}
					}
					else
					{ // children 4,5,6,7
						if (K->particles[index].x <= tree[j][k]->center.x)
						{ // children 4,7
							if (K->particles[index].y <= tree[j][k]->center.y)
							{ // child 4
								tree[J][Kp + 4]->chargeLocations.push_back(index);
							}
							else
							{ // child 7
								tree[J][Kp + 7]->chargeLocations.push_back(index);
							}
						}
						else
						{ // children 5,6
							if (K->particles[index].y <= tree[j][k]->center.y)
							{ // child 5
								tree[J][Kp + 5]->chargeLocations.push_back(index);
							}
							else
							{ // child 6
								tree[J][Kp + 6]->chargeLocations.push_back(index);
							}
						}
					}
				}
			}
		}
	}

	void print_charge_location1()
	{
		// for (int j = 0; j <= nLevels; j++) {
		for (int k = 0; k < nBoxesPerLevel[nLevels]; k++)
		{
			std::cout << tree[nLevels][k]->chargeLocations.size() << std::endl;
			// std::cout << "At " << "(" << j << " , " << k <<")" << std::endl;
			// for(auto& itr : tree[j][k]->chargeLocations){
			// 	std::cout << itr << "  ";
			// }
			// std::cout << std::endl;
			// std::cout << "=====================================================================" << std::endl;
		}
		std::cout << "=====================================================================" << std::endl;
		//}
	}

	void assignNonLeafChargeLocations()
	{
		for (int j = nLevels - 1; j >= 1; j--)
		{
			for (int k = 0; k < nBoxesPerLevel[j]; k++)
			{
				tree[j][k]->chargeLocations.clear();
				for (int c = 0; c < 8; c++)
				{
					tree[j][k]->chargeLocations.insert(tree[j][k]->chargeLocations.end(), tree[j + 1][8 * k + c]->chargeLocations.begin(), tree[j + 1][8 * k + c]->chargeLocations.end());
				}
			}
		}
	}

	void print_charge_location2()
	{
		for (int k = 0; k < nBoxesPerLevel[nLevels]; k++)
		{
			std::cout << tree[nLevels][k]->chargeLocations.size() << std::endl;
		}
	}

	void assignCharges(Vec &charges)
	{
		collective_charge = charges;
		for (int j = 1; j <= nLevels; j++)
		{
			for (int k = 0; k < nBoxesPerLevel[j]; k++)
			{
				tree[j][k]->charges = Vec::Zero(tree[j][k]->chargeLocations.size());
				for (size_t i = 0; i < tree[j][k]->chargeLocations.size(); i++)
				{
					tree[j][k]->charges[i] = charges[tree[j][k]->chargeLocations[i]];
				}
			}
		}
	}

	//! Far-field interaction using NCA

	void getNodes()
	{
		for (int j = nLevels; j >= 2; --j)
		{
			getNodes_incoming_level(j);
		}
	}

	void getParticlesFromChildren_incoming(int j, int k, std::vector<int> &searchNodes)
	{
		if (j == nLevels)
		{
			searchNodes.insert(searchNodes.end(), tree[j][k]->chargeLocations.begin(), tree[j][k]->chargeLocations.end());
		}
		else
		{
			int J = j + 1;
			for (int c = 0; c < 8; c++)
			{
				searchNodes.insert(searchNodes.end(), tree[J][8 * k + c]->incoming_checkPoints.begin(), tree[J][8 * k + c]->incoming_checkPoints.end());
			}
		}
	}

	void getNodes_incoming_box(int j, int k, int &ComputedRank)
	{
		std::vector<int> boxA_Nodes;
		getParticlesFromChildren_incoming(j, k, boxA_Nodes);
		std::vector<int> IL_Nodes; // indices
		for (size_t in = 0; in < tree[j][k]->interactionList.size(); ++in)
		{
			int kIL = tree[j][k]->interactionList[in];
			if (!tree[j][k]->isVertex[kIL])
			{
				std::vector<int> chargeLocations;
				getParticlesFromChildren_incoming(j, kIL, chargeLocations);
				IL_Nodes.insert(IL_Nodes.end(), chargeLocations.begin(), chargeLocations.end());
			}
		}
		std::vector<int> row_bases, col_bases;
		Mat Ac, Ar, L, R;
		LowRank *LR = new LowRank(K, TOL_POW, boxA_Nodes, IL_Nodes);
		LR->ACA_only_nodes(row_bases, col_bases, ComputedRank, Ac, Ar, L, R);
		if (ComputedRank > 0)
		{
			for (size_t r = 0; r < row_bases.size(); ++r)
			{
				tree[j][k]->incoming_checkPoints.push_back(boxA_Nodes[row_bases[r]]);
			}
			for (size_t c = 0; c < col_bases.size(); ++c)
			{
				tree[j][k]->incoming_chargePoints.push_back(IL_Nodes[col_bases[c]]);
			}

			Mat temp = R.triangularView<Eigen::Upper>().solve<Eigen::OnTheRight>(Ac);
			tree[j][k]->L2P[0] = L.triangularView<Eigen::Lower>().solve<Eigen::OnTheRight>(temp);
		}
		delete LR;
	}

	void getNodes_incoming_level(int j)
	{
		int ComputedRank;
		for (int k = 0; k < nBoxesPerLevel[j]; ++k)
		{
			getNodes_incoming_box(j, k, ComputedRank);
		}
	}

	void evaluate_M2M_far()
	{
		for (int j = nLevels; j >= 2; --j)
		{
			#pragma omp parallel for
			for (int k = 0; k < nBoxesPerLevel[j]; ++k)
			{
				Vec source_densities;
				if (j == nLevels)
				{
					source_densities = tree[j][k]->charges;
				}
				else
				{
					int J = j + 1;
					int Veclength = 0;
					for (int child = 0; child < 8; child++)
					{
						Veclength += tree[J][8 * k + child]->outgoing_charges.size();
					}
					source_densities = Vec::Zero(Veclength); // = tree[j][k]->multipoles//source densities
					int start = 0;
					for (int child = 0; child < 8; child++)
					{
						int NumElem = tree[J][8 * k + child]->outgoing_charges.size();
						source_densities.segment(start, NumElem) = tree[J][8 * k + child]->outgoing_charges;
						start += NumElem;
					}
				}
				clock_t start;
				start = clock();
				tree[j][k]->outgoing_charges = tree[j][k]->L2P[0].transpose() * source_densities;
				elapsed_mvp += (clock() - start) / (double)CLOCKS_PER_SEC;
			}
		}
	}

	void assemble_M2L_far(){
		for (int j = 2; j <= nLevels; ++j)
		{
			for (int k = 0; k < nBoxesPerLevel[j]; ++k)
			{
				for (size_t in = 0; in < tree[j][k]->interactionList.size(); ++in)
				{
					int kIL = tree[j][k]->interactionList[in];
					if (!tree[j][k]->isVertex[kIL])
					{
						if (tree[j][k]->M2L_far[kIL].size() == 0){
							tree[j][k]->M2L_far[kIL] = K->getMatrix(tree[j][k]->incoming_checkPoints, tree[j][kIL]->incoming_checkPoints);
						}
						if (tree[j][kIL]->M2L_far[k].size() == 0){
							tree[j][kIL]->M2L_far[k] = tree[j][k]->M2L_far[kIL].transpose();
						}
					}
				}
			}
		}
	}

	void evaluate_M2L_far()
	{
		for (int j = 2; j <= nLevels; ++j)
		{
			#pragma omp parallel for
			for (int k = 0; k < nBoxesPerLevel[j]; ++k)
			{
				tree[j][k]->incoming_potential = Vec::Zero(tree[j][k]->incoming_checkPoints.size());
				for (size_t in = 0; in < tree[j][k]->interactionList.size(); ++in)
				{
					int kIL = tree[j][k]->interactionList[in];
					if (!tree[j][k]->isVertex[kIL])
					{
						clock_t start;
						start = clock();
						tree[j][k]->incoming_potential += tree[j][k]->M2L_far[kIL] * tree[j][kIL]->outgoing_charges;
						elapsed_mvp += (clock() - start) / (double)CLOCKS_PER_SEC;
					}
				}
			}
		}
	}

	void evaluate_L2L_far()
	{
		for (int j = 2; j <= nLevels; ++j)
		{
			#pragma omp parallel for
			for (int k = 0; k < nBoxesPerLevel[j]; ++k)
			{
				if (j != nLevels)
				{
					clock_t start;
					start = clock();
					Vec temp = tree[j][k]->L2P[0] * tree[j][k]->incoming_potential;
					double duration = (clock() - start) / (double)CLOCKS_PER_SEC;
					elapsed_mvp += duration;

					int ind = 0;
					for (int c = 0; c < 8; ++c)
					{
						tree[j + 1][8 * k + c]->incoming_potential += temp.segment(ind, tree[j + 1][8 * k + c]->incoming_checkPoints.size());
						ind += tree[j + 1][8 * k + c]->incoming_checkPoints.size();
					}
				}
				else
				{
					clock_t start;
					start = clock();
					tree[j][k]->potential = tree[j][k]->L2P[0] * tree[j][k]->incoming_potential; // Local to particle
					double duration = (clock() - start) / (double)CLOCKS_PER_SEC;
					elapsed_mvp += duration;
				}
			}
		}
	}

	//! Vertex sharing interaction
	/////////////////////////////////////////////////////////////
	/*********** potential for vertex sharng interaction *************/
	////////////////////////////////////////////////////////////
	// Top down approach
	void getNodes_incoming_box_vertex(int j, int k, int &ComputedRank)
	{
		std::vector<int> boxA_Nodes = tree[j][k]->chargeLocations; // row
		std::vector<int> IL_Nodes;								  // col
		for (size_t i = 0; i < tree[j][k]->interactionList.size(); i++)
		{
			int ki = tree[j][k]->interactionList[i];
			if (tree[j][k]->isVertex[ki])
			{
				IL_Nodes.insert(IL_Nodes.end(), tree[j][ki]->chargeLocations.begin(), tree[j][ki]->chargeLocations.end());
			}
		}
		// Add the parent pivots except the root level
		if (j > 1)
		{
			IL_Nodes.insert(IL_Nodes.end(), tree[j - 1][k / 8]->incoming_chargePoints_ver.begin(), tree[j - 1][k / 8]->incoming_chargePoints_ver.end());
		}
		std::vector<int> row_bases, col_bases;
		Mat Ac, Ar;
		LowRank *LR = new LowRank(K, TOL_POW, boxA_Nodes, IL_Nodes);
		LR->ACA_only_nodes(row_bases, col_bases, ComputedRank, Ac, Ar, tree[j][k]->L_ver[0], tree[j][k]->R_ver[0]);
		if (ComputedRank > 0)
		{
			for (size_t r = 0; r < row_bases.size(); ++r)
			{
				tree[j][k]->incoming_checkPoints_ver.push_back(boxA_Nodes[row_bases[r]]);
			}
			for (size_t c = 0; c < col_bases.size(); ++c)
			{
				tree[j][k]->incoming_chargePoints_ver.push_back(IL_Nodes[col_bases[c]]);
			}
		}
		delete LR;
		// std::cout << "(" << j << "," << k << ")" << " Row = " << boxA_Nodes.size() << " Col = " << IL_Nodes.size() << " Rank = " << ComputedRank << std::endl;
	}

	void getNodes_ver()
	{
		for (int j = 1; j <= nLevels; ++j)
		{
			int ComputedRank;
			for (int k = 0; k < nBoxesPerLevel[j]; ++k)
			{
				getNodes_incoming_box_vertex(j, k, ComputedRank);
			}
		}
	}

	// Get M2M operator for vertex
	void assemble_M2M_ver()
	{
		for (int j = 1; j <= nLevels; ++j)
		{
			for (int k = 0; k < nBoxesPerLevel[j]; ++k)
			{
				if (j == nLevels)
				{
					Mat temp = K->getMatrix(tree[j][k]->incoming_chargePoints_ver, tree[j][k]->chargeLocations);
					Mat temp2 = tree[j][k]->R_ver[0].transpose().triangularView<Eigen::Lower>().solve(temp);
					tree[j][k]->M2M_ver[0] = tree[j][k]->L_ver[0].transpose().triangularView<Eigen::Upper>().solve(temp2);
				}
				else
				{
					for (int c = 0; c < 8; c++)
					{
						Mat temp = K->getMatrix(tree[j][k]->incoming_chargePoints_ver, tree[j + 1][8 * k + c]->incoming_checkPoints_ver);
						Mat temp2 = tree[j][k]->R_ver[0].transpose().triangularView<Eigen::Lower>().solve(temp);
						tree[j][k]->M2M_ver[c] = tree[j][k]->L_ver[0].transpose().triangularView<Eigen::Upper>().solve(temp2);
					}
				}
			}
		}
	}

	void evaluate_M2M_ver()
	{
		for (int j = nLevels; j >= 1; --j)
		{
			#pragma omp parallel for
			for (int k = 0; k < nBoxesPerLevel[j]; ++k)
			{
				if (j == nLevels)
				{
					clock_t start;
					start = clock();
					tree[j][k]->outgoing_charges = tree[j][k]->M2M_ver[0] * tree[j][k]->charges;
					elapsed_mvp += (clock() - start) / (double)CLOCKS_PER_SEC;
				}
				else
				{
					tree[j][k]->outgoing_charges = Vec::Zero(tree[j][k]->incoming_checkPoints_ver.size());
					for (int c = 0; c < 8; c++)
					{
						clock_t start;
						start = clock();
						tree[j][k]->outgoing_charges += tree[j][k]->M2M_ver[c] * tree[j + 1][8 * k + c]->outgoing_charges;
						elapsed_mvp += (clock() - start) / (double)CLOCKS_PER_SEC;
					}
				}
			}
		}
	}

	void assemble_M2L_ver(){
		for (int j = 1; j <= nLevels; ++j)
		{
			for (int k = 0; k < nBoxesPerLevel[j]; ++k)
			{
				tree[j][k]->incoming_potential_ver = Vec::Zero(tree[j][k]->incoming_checkPoints_ver.size());
				for (size_t i = 0; i < tree[j][k]->interactionList.size(); i++)
				{
					int kIL = tree[j][k]->interactionList[i];
					if (tree[j][k]->isVertex[kIL])
					{
						if (tree[j][k]->M2L_ver[kIL].size() == 0){
							tree[j][k]->M2L_ver[kIL] = K->getMatrix(tree[j][k]->incoming_checkPoints_ver, tree[j][kIL]->incoming_checkPoints_ver);
						}
						if (tree[j][kIL]->M2L_ver[k].size() == 0){
							tree[j][kIL]->M2L_ver[k] = tree[j][k]->M2L_ver[kIL].transpose();
						}
					}
				}
			}
		}
	}

	void evaluate_M2L_ver()
	{
		for (int j = 1; j <= nLevels; ++j)
		{
			#pragma omp parallel for
			for (int k = 0; k < nBoxesPerLevel[j]; ++k)
			{
				tree[j][k]->incoming_potential_ver = Vec::Zero(tree[j][k]->incoming_checkPoints_ver.size());
				for (size_t i = 0; i < tree[j][k]->interactionList.size(); i++)
				{
					int kIL = tree[j][k]->interactionList[i];
					if (tree[j][k]->isVertex[kIL])
					{
						clock_t start;
						start = clock();
						tree[j][k]->incoming_potential_ver += tree[j][k]->M2L_ver[kIL] * tree[j][kIL]->outgoing_charges;
						elapsed_mvp += (clock() - start) / (double)CLOCKS_PER_SEC;
					}
				}
			}
		}
	}

	void evaluate_L2L_ver()
	{
		for (int j = 1; j <= nLevels; ++j)
		{
			#pragma omp parallel for
			for (int k = 0; k < nBoxesPerLevel[j]; ++k)
			{
				if (j != nLevels)
				{
					for (int c = 0; c < 8; c++)
					{

						clock_t start;
						start = clock();
						tree[j + 1][8 * k + c]->incoming_potential_ver += tree[j][k]->M2M_ver[c].transpose() * tree[j][k]->incoming_potential_ver;
						elapsed_mvp += (clock() - start) / (double)CLOCKS_PER_SEC;
					}
				}
				else
				{
					clock_t start;
					start = clock();
					tree[j][k]->potential_ver = tree[j][k]->M2M_ver[0].transpose() * tree[j][k]->incoming_potential_ver;
					elapsed_mvp += (clock() - start) / (double)CLOCKS_PER_SEC;
				}
			}
		}
	}

    void getUVtTree_nn() {
        for (int j = 1; j <= nLevels; j++) {
            for (int k = 0; k < nBoxesPerLevel[j]; k++) {
                for (size_t i = 0; i < tree[j][k]->interactionList.size(); i++) {
                    int ki = tree[j][k]->interactionList[i];
                    getUVtInstance(j,k,ki,tree[j][k]->L[i],tree[j][k]->R[i], tree[j][k]->row_basis[i],tree[j][k]->col_basis[i]);
                }
            }
        }
    }

	void getUVtTree_sn()
	{
		for (int j = 1; j <= nLevels; j++)
		{
			for (int k = 0; k < nBoxesPerLevel[j]; k++)
			{
				for (size_t i = 0; i < tree[j][k]->interactionList.size(); i++)
				{
					int ki = tree[j][k]->interactionList[i];
					if (tree[j][k]->isVertex[ki])
					{
						getUVtInstance(j, k, ki, tree[j][k]->L[i], tree[j][k]->R[i], tree[j][k]->row_basis[i], tree[j][k]->col_basis[i]);
					}
				}
			}
		}
	}

    void getUVtInstance(int j, int k, int ki, Mat& L,Mat& R, std::vector<int>& row_indices, std::vector<int>& col_indices) {
        int computed_rank;
        std::vector<int> row_indices_local,col_indices_local;
        LowRank* LR             =       new LowRank(K, TOL_POW, tree[j][k]->chargeLocations, tree[j][ki]->chargeLocations);
        LR->ACA_only_nodesCUR(row_indices_local, col_indices_local, computed_rank, L, R);
        for (int r = 0; r < computed_rank; r++) {
            row_indices.push_back(tree[j][k]->chargeLocations[row_indices_local[r]]);
        }
        for (int c = 0; c < computed_rank; c++) {
            col_indices.push_back(tree[j][ki]->chargeLocations[col_indices_local[c]]);
        }
        delete LR;
    }

    void assembleAcAr_nn(){
        for (int j = 1; j <= nLevels; j++) {
            for (int k = 0; k < nBoxesPerLevel[j]; k++) {
                for (size_t i = 0; i < tree[j][k]->interactionList.size(); i++) {
                    if (tree[j][k]->col_basis[i].size() > 0) {
                        int ki = tree[j][k]->interactionList[i];
                        if(tree[j][k]->Ac[ki].size() == 0)
                            tree[j][k]->Ac[ki] = K->getMatrix(tree[j][k]->chargeLocations,tree[j][k]->col_basis[i]);
                        if(tree[j][k]->Ar[ki].size() == 0)
                        tree[j][k]->Ar[ki] = K->getMatrix(tree[j][k]->row_basis[i],tree[j][ki]->chargeLocations);
                    }
            	}
            }
        }
    }

	void assembleAcAr_sn()
	{
		for (int j = 1; j <= nLevels; j++)
		{
			for (int k = 0; k < nBoxesPerLevel[j]; k++)
			{
				for (size_t i = 0; i < tree[j][k]->interactionList.size(); i++)
				{
					if (tree[j][k]->col_basis[i].size() > 0)
					{
						int ki = tree[j][k]->interactionList[i];

						if (tree[j][k]->isVertex[ki])
						{
							if (tree[j][k]->Ac[ki].size() == 0)
								tree[j][k]->Ac[ki] = K->getMatrix(tree[j][k]->chargeLocations, tree[j][k]->col_basis[i]);
							if (tree[j][k]->Ar[ki].size() == 0)
								tree[j][k]->Ar[ki] = K->getMatrix(tree[j][k]->row_basis[i], tree[j][ki]->chargeLocations);
						}
					}
				}
			}
		}
	}

	void evaluateFarField_sn()
	{
		for (int j = 1; j <= nLevels; j++)
		{
			#pragma omp parallel for
			for (int k = 0; k < nBoxesPerLevel[j]; k++)
			{
				tree[j][k]->potential_ver = Vec::Zero(tree[j][k]->charges.size());
				for (size_t i = 0; i < tree[j][k]->interactionList.size(); i++)
				{
					if (tree[j][k]->col_basis[i].size() > 0)
					{
						int ki = tree[j][k]->interactionList[i];

						if (tree[j][k]->isVertex[ki])
						{
							clock_t start5;
							start5 = clock();
							Vec t0 = tree[j][k]->Ar[ki] * tree[j][ki]->charges;
							Vec t1 = tree[j][k]->L[i].triangularView<Eigen::Lower>().solve(t0);
							Vec t2 = tree[j][k]->R[i].triangularView<Eigen::Upper>().solve(t1);
							tree[j][k]->potential_ver += tree[j][k]->Ac[ki] * t2;
							elapsed_mvp += (clock() - start5) / (double)CLOCKS_PER_SEC;
						}
					}
				}
			}
		}
	}

    void evaluateFarField_nn() {
        for (int j = 1; j <= nLevels; j++) {
			#pragma omp parallel for
            for (int k = 0; k < nBoxesPerLevel[j]; k++) {
                tree[j][k]->potential = Vec::Zero(tree[j][k]->charges.size());
                for (size_t i = 0; i < tree[j][k]->interactionList.size(); i++) {
                    if (tree[j][k]->col_basis[i].size() > 0) {
                        int ki = tree[j][k]->interactionList[i];
						
						clock_t start5;
						start5 = clock();
                        Vec t0 = tree[j][k]->Ar[ki]*tree[j][ki]->charges;
                        Vec t1 = tree[j][k]->L[i].triangularView<Eigen::Lower>().solve(t0);
                        Vec t2 = tree[j][k]->R[i].triangularView<Eigen::Upper>().solve(t1);
                        tree[j][k]->potential += tree[j][k]->Ac[ki]*t2;
						elapsed_mvp += (clock() - start5) / (double)CLOCKS_PER_SEC;
                    }
            	}
            }
        }
    }

	void assemble_NearField()
	{ // evaluating at chargeLocations
		#pragma omp parallel for
		for (int k = 0; k < nBoxesPerLevel[nLevels]; k++)
		{
			for (int n = 0; n < 26; n++)
			{
				int nn = tree[nLevels][k]->neighborNumbers[n];
				if (nn != -1)
				{
					tree[nLevels][k]->denseMatrices[n] = K->getMatrix(tree[nLevels][k]->chargeLocations, tree[nLevels][nn]->chargeLocations);;
				}
			}
			tree[nLevels][k]->denseMatrices[26] = K->getMatrix(tree[nLevels][k]->chargeLocations, tree[nLevels][k]->chargeLocations);;
		}
	}

	void evaluate_NearField()
	{
        #if defined(USE_nHODLRdD) || defined(USE_snHODLRdD)
        if (nLevels < 2)
        {
            for (int k = 0; k < nBoxesPerLevel[nLevels]; ++k)
            {
                tree[nLevels][k]->potential = Vec::Zero(tree[nLevels][k]->charges.size());
            }
        }
        #endif
		#pragma omp parallel for
		for (int k = 0; k < nBoxesPerLevel[nLevels]; k++)
		{
			for (int n = 0; n < 26; n++)
			{
				int nn = tree[nLevels][k]->neighborNumbers[n];
				if (nn != -1)
				{
					clock_t start;
					start = clock();
					tree[nLevels][k]->potential += tree[nLevels][k]->denseMatrices[n]  *  tree[nLevels][nn]->charges;
					elapsed_mvp += (clock() - start) / (double)CLOCKS_PER_SEC;
				}
			}

			clock_t start;
			start = clock();
			tree[nLevels][k]->potential += tree[nLevels][k]->denseMatrices[26] * tree[nLevels][k]->charges; // self Interaction
			elapsed_mvp += (clock() - start) / (double)CLOCKS_PER_SEC;
		}
	}

	// The computed total potential
	void collect_all_potentials_n(Vec &potential)
	{
		potential = Vec::Zero(N);
		for (int j = 1; j <= nLevels; j++)
		{
			int count = 0;
			Vec potentialTemp = Vec::Zero(N);
			for (int k = 0; k < nBoxesPerLevel[j]; k++)
			{ // using the fact that all the leaves have same number of particles
				potentialTemp.segment(count, tree[j][k]->potential.size()) = tree[j][k]->potential + tree[j][k]->potential_ver;
				count += tree[j][k]->potential.size();
			}
			potential += potentialTemp;
		}
	}

	void collect_all_potentials_semi_1(Vec &potential)
	{
		potential = Vec::Zero(N);
		for (int j = 1; j <= nLevels; j++)
		{
			int count = 0;
			Vec potentialTemp = Vec::Zero(N);
			for (int k = 0; k < nBoxesPerLevel[j]; k++)
			{
				potentialTemp.segment(count, tree[j][k]->potential.size()) = tree[j][k]->potential;
				count += tree[j][k]->potential.size();
			}
			potential += potentialTemp;
		}
	}

	void collect_all_potentials_semi_2(Vec &potential)
	{
		for (int j = 1; j <= nLevels; j++)
		{
			int count = 0;
			Vec potentialTemp = Vec::Zero(N);
			for (int k = 0; k < nBoxesPerLevel[j]; k++)
			{
				potentialTemp.segment(count, tree[j][k]->potential_ver.size()) = tree[j][k]->potential_ver;
				count += tree[j][k]->potential_ver.size();
			}
			potential += potentialTemp;
		}
	}

	void collect_all_potentials_nn(Vec &potential)
	{
		potential = Vec::Zero(N);
		for (int j = 1; j <= nLevels; j++)
		{
			int count = 0;
			Vec potentialTemp = Vec::Zero(N);
			for (int k = 0; k < nBoxesPerLevel[j]; k++)
			{ // using the fact that all the leaves have same number of particles
				potentialTemp.segment(count, tree[j][k]->potential.size()) = tree[j][k]->potential;
				count += tree[j][k]->potential.size();
			}
			potential += potentialTemp;
		}
	}

	void reorder(Vec &potential)
	{
		Vec potentialTemp = potential;
		int count = 0;
		// std::cout << "index: " << std::endl;
		for (int k = 0; k < nBoxesPerLevel[nLevels]; ++k)
		{
			for (size_t i = 0; i < tree[nLevels][k]->chargeLocations.size(); ++i)
			{
				int index = tree[nLevels][k]->chargeLocations[i];
				// std::cout << index << std::endl;
				potential(index) = potentialTemp(count);
				++count;
			}
		}
		collective_potential = potential;
	}

	void findMemory_n() {
		for (int j = 2; j <= nLevels; j++) {
			for (int k = 0; k < nBoxesPerLevel[j]; k++) {
				memory += 2*tree[j][k]->L2P[0].size();
			}
		}

		for (int j = 1; j <= nLevels; j++) {
			for (int k = 0; k < nBoxesPerLevel[j]; k++) {
				if (j == nLevels)
				{
					memory += 2*tree[j][k]->M2M_ver[0].size();
				}
				else
				{
					for (int c = 0; c < 8; c++)
					{
						memory += 2*tree[j][k]->M2M_ver[c].size();
					}
				}
			}
		}

		// M2L memory
		for (int j = 2; j <= nLevels; j++) {
			for (int k = 0; k < nBoxesPerLevel[j]; k++) {
				for (size_t i = 0; i < tree[j][k]->interactionList.size(); i++) {
					int kIL = tree[j][k]->interactionList[i];
					if (!tree[j][k]->isVertex[kIL]){
						memory += tree[j][k]->M2L_far[kIL].size();
					}
				}
			}
		}
		for (int j = 1; j <= nLevels; j++) {
			for (int k = 0; k < nBoxesPerLevel[j]; k++) {
				for (size_t i = 0; i < tree[j][k]->interactionList.size(); i++) {
					int kIL = tree[j][k]->interactionList[i];
					if (tree[j][k]->isVertex[kIL]){
						memory += tree[j][k]->M2L_ver[kIL].size();
					}
				}
			}
		}
		for (int k = 0; k < nBoxesPerLevel[nLevels]; k++) {
			memory += tree[nLevels][k]->chargeLocations.size()*tree[nLevels][k]->chargeLocations.size(); //self
			for (size_t n = 0; n < 26; n++) {
				int nn = tree[nLevels][k]->neighborNumbers[n];
				if(nn != -1) {
					memory += tree[nLevels][k]->chargeLocations.size()*tree[nLevels][nn]->chargeLocations.size();
				}
			}
		}
	}

	void findMemory_nn() {
		for (int j = 1; j <= nLevels; j++) {
			for (int k = 0; k < nBoxesPerLevel[j]; k++) {
				for (size_t i = 0; i < tree[j][k]->interactionList.size(); i++) {
					int ki = tree[j][k]->interactionList[i];
					memory += tree[j][k]->L[i].size();
					memory += tree[j][k]->R[i].size();
					memory += tree[j][k]->Ac[ki].size();
					memory += tree[j][k]->Ar[ki].size();
				}
			}
		}
		for (int k = 0; k < nBoxesPerLevel[nLevels]; k++) {
			memory += tree[nLevels][k]->chargeLocations.size()*tree[nLevels][k]->chargeLocations.size(); //self
			for (int n = 0; n < 26; n++) {
				int nn = tree[nLevels][k]->neighborNumbers[n];
				if(nn != -1) {
					memory += tree[nLevels][k]->chargeLocations.size()*tree[nLevels][nn]->chargeLocations.size();
				}
			}
		}
	}

	void findMemory_sn() {
		for (int j = 2; j <= nLevels; j++) {
			for (int k = 0; k < nBoxesPerLevel[j]; k++) {
				memory += 2*tree[j][k]->L2P[0].size();
			}
		}

		// M2L memory
		for (int j = 2; j <= nLevels; j++) {
			for (int k = 0; k < nBoxesPerLevel[j]; k++) {
				for (size_t i = 0; i < tree[j][k]->interactionList.size(); i++) {
					int kIL = tree[j][k]->interactionList[i];
					if (!tree[j][k]->isVertex[kIL]){
						memory += tree[j][k]->M2L_far[kIL].size();
					}
				}
			}
		}
		for (int k = 0; k < nBoxesPerLevel[nLevels]; k++) {
			memory += tree[nLevels][k]->chargeLocations.size()*tree[nLevels][k]->chargeLocations.size(); //self
			for (size_t n = 0; n < 26; n++) {
				int nn = tree[nLevels][k]->neighborNumbers[n];
				if(nn != -1) {
					memory += tree[nLevels][k]->chargeLocations.size()*tree[nLevels][nn]->chargeLocations.size();
				}
			}
		}
		for (int j = 1; j <= nLevels; j++) {
			for (int k = 0; k < nBoxesPerLevel[j]; k++) {
				for (size_t i = 0; i < tree[j][k]->interactionList.size(); i++) {
					int ki = tree[j][k]->interactionList[i];
					if (tree[j][k]->isVertex[ki]){
						memory += tree[j][k]->L[i].size();
						memory += tree[j][k]->R[i].size();
						memory += tree[j][k]->Ac[ki].size();
						memory += tree[j][k]->Ar[ki].size();
					}
				}
			}
		}
	}

	Vec getTruePoten(Vec& charge){
		Vec true_poten(N);
		std::vector<int> v;
		for (int i = 0; i < N; ++i)
		{
			v.push_back(i);
		}
		#pragma omp parallel for
		for (int i = 0; i < N; ++i)
		{
			true_poten(i) = K->getRow(i, v).transpose() * charge;
		}
		return true_poten;
	}

	// Total error
	double compute_total_error()
	{
		Vec true_poten(N);
		std::vector<int> v;
		for (int i = 0; i < N; ++i)
		{
			v.push_back(i);
		}
		#pragma omp parallel for
		for (int i = 0; i < N; ++i)
		{
			true_poten(i) = K->getRow(i, v).transpose() * collective_charge;
		}
		Vec errVec = true_poten - collective_potential;
		double error = errVec.norm() / true_poten.norm();
		return error;
	}
};

#endif
