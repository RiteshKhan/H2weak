
#ifndef TREE_HPP
#define TREE_HPP

#include "kernel.hpp"
#include "HODLR1Dnode.hpp"
#include "aca.hpp"

template <typename kerneltype>
class hodlr1d_tree
{
public:
    kerneltype *K;
    int nLevels; //	Number of levels in the tree.
    int N;       //	Number of particles.
    double L;    //	Semi-length of the simulation box.
    int nParticlesInLeaf;
    int TOL_POW;

    std::vector<int> nNodesPerLevel;             //	Number of boxes per level in the tree
    std::vector<double> nodeLength;              //	Line radius per level in the tree assuming the line at the root is [-1,1]
    std::vector<std::vector<hodlr_node *>> tree; //	The tree storing all the information.
    std::vector<double> Nodes;
    std::vector<double> gridPoints; // all particles in domain
    Vec collective_charge;
    Vec collective_potential, collective_potential_nn;

    double elapsed_mvp = 0.0;
    double elapsed_assem = 0.0;
    double memory = 0.0;

    hodlr1d_tree(kerneltype *K, int N, int nLevels, int nParticlesInLeaf, double L, int TOL_POW)
    {
        this->K = K;
        this->nLevels = nLevels;
        this->L = L;
        this->nParticlesInLeaf = nParticlesInLeaf;
        this->TOL_POW = TOL_POW;
        this->N = N;
        nNodesPerLevel.push_back(1);
        nodeLength.push_back(L);
        for (int k = 1; k <= nLevels; ++k)
        {
            nNodesPerLevel.push_back(2 * nNodesPerLevel[k - 1]);
            nodeLength.push_back(0.5 * nodeLength[k - 1]);
        }
    }

    // Chebyshev nodes
    void set_Standard_Cheb_Nodes()
    {
        for (int k = 0; k < N; ++k)
        {
            K->particles.push_back(-cos((k + 0.5) / N * PI));
        }
    }

    // Uniform nodes
    void set_Uniform_Nodes()
    {
        for (int k = 0; k < N; ++k)
        {
            K->particles.push_back(-1.0 + 2.0 * (k + 1.0) / (N + 1.0));
        }
    }

    // Need to shift location at leaf level
    void shift_Nodes(double radius, double center, std::vector<double> &particle_loc)
    {
        for (size_t i = 0; i < Nodes.size(); i++)
        {
            double temp;
            temp = Nodes[i] * radius + center;
            particle_loc.push_back(temp);
        }
    }

    // Create tree and interaction list
    void createTree()
    {
        //	First create root and add to tree
        hodlr_node *root = new hodlr_node;
        root->nodeNumber = 0;
        root->parentNumber = -1;
        for (int l = 0; l < 2; ++l)
        {
            root->childrenNumbers[l] = l;
        }
        for (int l = 0; l < 2; ++l)
        {
            root->neighborNumbers[l] = -1;
        }
        std::vector<hodlr_node *> rootLevel;
        rootLevel.push_back(root);
        tree.push_back(rootLevel);

        for (int j = 1; j <= nLevels; ++j)
        {
            std::vector<hodlr_node *> level;
            for (int k = 0; k < nNodesPerLevel[j]; ++k)
            {
                hodlr_node *temp = new hodlr_node;
                temp->nodeNumber = k;
                temp->parentNumber = k / 2;
                for (int l = 0; l < 2; ++l)
                {
                    temp->childrenNumbers[l] = 2 * k + l;
                }
                level.push_back(temp);
            }
            tree.push_back(level);
        }
    }

    // Assigns the interactions for child 0 of a node
    void assign_Child0_Interaction(int j, int k)
    {
        int nL = j + 1;
        int nC = 2 * k;
        // Assigns siblings
        // ____**____|____N1____
        tree[nL][nC]->neighborNumbers[1] = nC + 1;
    }

    // Assigns the interactions for child 1 of a line
    void assign_Child1_Interaction(int j, int k)
    {
        int nL = j + 1;
        int nC = 2 * k + 1;
        // Assigns siblings
        // ____N0____|____**____
        tree[nL][nC]->neighborNumbers[0] = nC - 1;
    }

    // Interaction of all levels accross the tree
    void assign_all_interaction()
    {
        for (int j = 0; j < nLevels; ++j)
        {
            for (int k = 0; k < nNodesPerLevel[j]; ++k)
            {
                assign_Child0_Interaction(j, k);
                assign_Child1_Interaction(j, k);
            }
        }
    }

    // Center location of different nodes
    void assign_center_location()
    {
        int J;
        tree[0][0]->center = 0.0;
        for (int j = 0; j < nLevels; ++j)
        {
            J = j + 1;
            double shift = 0.5 * nodeLength[j];
            #pragma omp parallel for
            for (int k = 0; k < nNodesPerLevel[j]; ++k)
            {
                tree[J][2 * k]->center = tree[j][k]->center - shift;
                tree[J][2 * k + 1]->center = tree[j][k]->center + shift;
            }
        }
    }

    // Assign charge location
    void assignChargeLocations()
    {
        for (int i = 0; i < N; i++)
        {
            tree[0][0]->chargeLocations.push_back(i);
        }
        for (int j = 0; j < nLevels; j++)
        { // assign particles to its children
            for (int k = 0; k < nNodesPerLevel[j]; k++)
            {
                int J = j + 1;
                int Kp = 2 * k;
                for (size_t i = 0; i < tree[j][k]->chargeLocations.size(); i++)
                {
                    int index = tree[j][k]->chargeLocations[i];
                    if (K->particles[index] <= tree[j][k]->center)
                    {
                        // child 0
                        tree[J][Kp]->chargeLocations.push_back(index);
                    }
                    else
                    {
                        // child 1
                        tree[J][Kp + 1]->chargeLocations.push_back(index);
                    }
                }
            }
        }
    }

    // Assign charge at non leaf
    void assignNonLeafChargeLocations()
    {
        for (int j = nLevels - 1; j > 0; --j)
        {
            for (int k = 0; k < nNodesPerLevel[j]; ++k)
            {
                tree[j][k]->chargeLocations.clear();
                for (int c = 0; c < 2; ++c)
                {
                    tree[j][k]->chargeLocations.insert(tree[j][k]->chargeLocations.end(), tree[j + 1][2 * k + c]->chargeLocations.begin(), tree[j + 1][2 * k + c]->chargeLocations.end());
                }
            }
        }
    }

    // Assign charges
    void assignCharges(Vec &charges)
    {
        collective_charge = charges;
        for (int j = 1; j <= nLevels; j++)
        {
            for (int k = 0; k < nNodesPerLevel[j]; k++)
            {
                tree[j][k]->charges = Vec::Zero(tree[j][k]->chargeLocations.size());
                for (size_t i = 0; i < tree[j][k]->chargeLocations.size(); i++)
                {
                    tree[j][k]->charges[i] = charges[tree[j][k]->chargeLocations[i]];
                }
            }
        }
    }

    void getNodes_incoming_box(int j, int k, int &ComputedRank)
    {
        std::vector<int> boxA_Nodes = tree[j][k]->chargeLocations; // row
        std::vector<int> IL_Nodes;                                 // col
        for (int nn = 0; nn < 2; ++nn)
        {
            int p = tree[j][k]->neighborNumbers[nn];
            if (p != -1)
            {
                IL_Nodes.insert(IL_Nodes.end(), tree[j][p]->chargeLocations.begin(), tree[j][p]->chargeLocations.end());
            }
        }
        // Add the parent pivots except the root level
        if (j > 1)
        {
            IL_Nodes.insert(IL_Nodes.end(), tree[j - 1][k / 2]->incoming_chargePoints.begin(), tree[j - 1][k / 2]->incoming_chargePoints.end());
        }
        std::vector<int> row_bases, col_bases;
        Mat Ac, Ar;
        LowRank *LR = new LowRank(K, TOL_POW, boxA_Nodes, IL_Nodes);
        LR->ACA_only_nodes(row_bases, col_bases, ComputedRank, Ac, Ar, tree[j][k]->L, tree[j][k]->R);
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
        }
        // std::cout << "(" << j << "," << k << ")" << " Row = " << boxA_Nodes.size() << " Col = " << IL_Nodes.size() << " Rank = " << ComputedRank << std::endl;
    }

    void getNodes()
    {
        for (int j = 1; j <= nLevels; ++j)
        {
            int ComputedRank;
            for (int k = 0; k < nNodesPerLevel[j]; ++k)
            {
                getNodes_incoming_box(j, k, ComputedRank);
            }
        }
    }

    // Get M2M operator for vertex
    void assemble_M2M()
    {
        for (int j = 1; j <= nLevels; j++)
        {
            for (int k = 0; k < nNodesPerLevel[j]; k++)
            {
                if (j == nLevels)
                {
                    Mat temp = K->getMatrix(tree[j][k]->incoming_chargePoints, tree[j][k]->chargeLocations);
                    Mat temp2 = tree[j][k]->R.transpose().triangularView<Eigen::Lower>().solve(temp);
                    tree[j][k]->M2M[0] = tree[j][k]->L.transpose().triangularView<Eigen::Upper>().solve(temp2);
                }
                else
                {
                    for (int c = 0; c < 2; c++)
                    {
                        Mat temp = K->getMatrix(tree[j][k]->incoming_chargePoints, tree[j + 1][2 * k + c]->incoming_checkPoints);
                        Mat temp2 = tree[j][k]->R.transpose().triangularView<Eigen::Lower>().solve(temp);
                        tree[j][k]->M2M[c] = tree[j][k]->L.transpose().triangularView<Eigen::Upper>().solve(temp2);
                    }
                }
            }
        }
    }

    void evaluate_M2M()
    {
        for (int j = nLevels; j >= 1; j--)
        {
            #pragma omp parallel for
            for (int k = 0; k < nNodesPerLevel[j]; k++)
            {
                if (j == nLevels)
                {
                    clock_t start8;
                    start8 = clock();
                    Mat M2M = tree[j][k]->M2M[0];
                    Vec w = tree[j][k]->charges;
                    elapsed_assem += (clock() - start8) / (double)CLOCKS_PER_SEC;

                    clock_t start;
                    start = clock();
                    tree[j][k]->outgoing_charges = M2M * w;
                    elapsed_mvp += (clock() - start) / (double)CLOCKS_PER_SEC;
                }
                else
                {
                    tree[j][k]->outgoing_charges = Vec::Zero(tree[j][k]->incoming_checkPoints.size());
                    for (int c = 0; c < 2; c++)
                    {
                        clock_t start9;
                        start9 = clock();
                        Mat M2M = tree[j][k]->M2M[c];
                        Vec w = tree[j + 1][2 * k + c]->outgoing_charges;
                        elapsed_assem += (clock() - start9) / (double)CLOCKS_PER_SEC;

                        clock_t start;
                        start = clock();
                        tree[j][k]->outgoing_charges += M2M * w;
                        elapsed_mvp += (clock() - start) / (double)CLOCKS_PER_SEC;
                    }
                }
            }
        }
    }

    // Get M2L operator for vertex
    void assemble_M2L()
    {
        for (int j = 1; j <= nLevels; j++)
        {
            for (int k = 0; k < nNodesPerLevel[j]; k++)
            {
                for (int nn = 0; nn < 2; ++nn)
                {
                    int kIL = tree[j][k]->neighborNumbers[nn];
                    if (kIL != -1)
                    {
                        if (tree[j][k]->M2L[kIL].size() == 0)
                            tree[j][k]->M2L[kIL] = K->getMatrix(tree[j][k]->incoming_checkPoints, tree[j][kIL]->incoming_checkPoints);
                        if (tree[j][kIL]->M2L[k].size() == 0)
                            tree[j][kIL]->M2L[k] = tree[j][k]->M2L[kIL].transpose();
                    }
                }
            }
        }
    }

    void evaluate_M2L()
    {
        for (int j = 1; j <= nLevels; ++j)
        {
            #pragma omp parallel for
            for (int k = 0; k < nNodesPerLevel[j]; ++k)
            {
                tree[j][k]->incoming_potential = Vec::Zero(tree[j][k]->incoming_checkPoints.size());
                for (int nn = 0; nn < 2; ++nn)
                {
                    int kIL = tree[j][k]->neighborNumbers[nn];
                    if (kIL != -1)
                    {
                        clock_t start;
                        start = clock();
                        tree[j][k]->incoming_potential += tree[j][k]->M2L[kIL] * tree[j][kIL]->outgoing_charges;
                        elapsed_mvp += (clock() - start) / (double)CLOCKS_PER_SEC;
                    }
                }
            }
        }
    }

    void evaluate_L2L()
    {
        for (int j = 1; j <= nLevels; j++)
        {
            #pragma omp parallel for
            for (int k = 0; k < nNodesPerLevel[j]; k++)
            {
                if (j != nLevels)
                {
                    for (int c = 0; c < 2; c++)
                    {
                        clock_t start;
                        start = clock();
                        tree[j + 1][2 * k + c]->incoming_potential += tree[j][k]->M2M[c].transpose() * tree[j][k]->incoming_potential;
                        elapsed_mvp += (clock() - start) / (double)CLOCKS_PER_SEC;
                    }
                }
                else
                {
                    clock_t start;
                    start = clock();
                    tree[j][k]->potential = tree[j][k]->M2M[0].transpose() * tree[j][k]->incoming_potential;
                    elapsed_mvp += (clock() - start) / (double)CLOCKS_PER_SEC;
                }
            }
        }
    }

    // Assemble Near-field
    void assemble_NearField()
    {
        #pragma omp parallel for
        for (int k = 0; k < nNodesPerLevel[nLevels]; k++)
        {
            tree[nLevels][k]->denseMatrices = K->getMatrix(tree[nLevels][k]->chargeLocations, tree[nLevels][k]->chargeLocations);
        }
    }

    // Evaluate Near-field
    void evaluate_NearField()
    {
        #pragma omp parallel for
        for (int k = 0; k < nNodesPerLevel[nLevels]; ++k)
        {
            // Self interaction
            tree[nLevels][k]->potential += tree[nLevels][k]->denseMatrices * tree[nLevels][k]->charges;
        }
    }

    // The computed total potential
    void collect_all_potentials(Vec &potential)
    {
        potential = Vec::Zero(N);
        for (int j = 1; j <= nLevels; j++)
        {
            int count = 0;
            Vec potentialTemp = Vec::Zero(N);
            for (int k = 0; k < nNodesPerLevel[j]; k++)
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
        for (int k = 0; k < nNodesPerLevel[nLevels]; ++k)
        {
            for (size_t i = 0; i < tree[nLevels][k]->chargeLocations.size(); ++i)
            {
                int index = tree[nLevels][k]->chargeLocations[i];
                potential(index) = potentialTemp(count);
                ++count;
            }
        }
        collective_potential = potential;
    }

    // void findMemory() {
    // 	for (int j = 2; j <= nLevels; j++) {
    // 		for (int k = 0; k < nBoxesPerLevel[j]; k++) {
    // 			memory += 2*tree[j][k].L2P[0].size();
    // 		}
    // 	}

    // 	for (int j = 1; j <= nLevels; j++) {
    // 		for (int k = 0; k < nBoxesPerLevel[j]; k++) {
    // 			if (j == nLevels)
    // 			{
    // 				memory += 2*tree[j][k].M2M_ver[0].size();
    // 			}
    // 			else
    // 			{
    // 				for (int c = 0; c < 4; c++)
    // 				{
    // 					memory += 2*tree[j][k].M2M_ver[c].size();
    // 				}
    // 			}
    // 		}
    // 	}

    // 	// M2L memory
    // 	for (int j = 2; j <= nLevels; j++) {
    // 		for (int k = 0; k < nBoxesPerLevel[j]; k++) {
    //             for (int in = 0; in < 12; ++in)
    //             {
    //                 int kIL = tree[j][k].innerNumbers[in];
    //                 if (kIL != -1)
    //                 {
    //                     memory += tree[j][k].M2L_far[kIL].size();
    //                 }
    //             }
    //             for (int on = 0; on < 12; ++on)
    //             {
    //                 int kIL = tree[j][k].outerNumbers[on];
    //                 if (kIL != -1)
    //                 {
    //                     memory += tree[j][k].M2L_far[kIL].size();
    //                 }
    //             }
    // 		}
    // 	}
    // 	for (int j = 1; j <= nLevels; j++) {
    // 		for (int k = 0; k < nBoxesPerLevel[j]; k++) {
    //             for (int co = 0; co < 4; ++co)
    //             {
    //                 int kIL = tree[j][k].cornerNumbers[co];
    //                 if (kIL != -1)
    //                 {
    //                     memory += tree[j][k].M2L_ver[kIL].size();
    //                 }
    //             }
    // 		}
    // 	}
    // 	for (int k = 0; k < nBoxesPerLevel[nLevels]; k++) {
    // 		memory += tree[nLevels][k].chargeLocations.size()*tree[nLevels][k].chargeLocations.size(); //self
    // 		for (size_t n = 0; n < 4; n++) {
    // 			int nn = tree[nLevels][k].neighborNumbers[n];
    // 			if(nn != -1) {
    // 				memory += tree[nLevels][k].chargeLocations.size()*tree[nLevels][nn].chargeLocations.size();
    // 			}
    // 		}
    // 	}
    // }



    void getUVtInstance(int j, int k, int kIL, Mat &U, Mat &V)
    {
        int computed_rank;
        std::vector<int> row_bases, col_bases;
        LowRank *LR = new LowRank(K, TOL_POW, tree[j][k]->chargeLocations, tree[j][kIL]->chargeLocations);
        LR->ACA_only_nodesUV(row_bases, col_bases, computed_rank, tree[j][k]->U[kIL], tree[j][k]->V[kIL]);
        delete LR;
    }

    void getUVtTree_nn()
    {
        for (int j = 1; j <= nLevels; j++)
        {
            for (int k = 0; k < nNodesPerLevel[j]; k++)
            {
                for (int nn = 0; nn < 2; ++nn)
                {
                    int kIL = tree[j][k]->neighborNumbers[nn];
                    if (kIL != -1)
                    {
                        getUVtInstance(j, k, kIL, tree[j][k]->U[kIL], tree[j][k]->V[kIL]);
                    }
                }
            }
        }
    }

    void fast_multiplication_nn()
    {
        for (int j = 1; j <= nLevels; ++j)
        {
            #pragma omp parallel for
            for (int k = 0; k < nNodesPerLevel[j]; ++k)
            {
                tree[j][k]->potential = Vec::Zero(tree[j][k]->charges.size());
                for (int nn = 0; nn < 2; ++nn)
                {
                    int kIL = tree[j][k]->neighborNumbers[nn];
                    if (kIL != -1)
                    {
                        tree[j][k]->potential += tree[j][k]->U[kIL] * (tree[j][k]->V[kIL] * tree[j][kIL]->charges);
                    }
                }
            }
        }
    }

    void collect_all_potential_nn(Vec &potential)
    {
        potential = Vec::Zero(N);
        for (int j = 1; j <= nLevels; j++)
        {
            int count = 0;
            Vec potentialTemp = Vec::Zero(N);
            for (int k = 0; k < nNodesPerLevel[j]; k++)
            { // using the fact that all the leaves have same number of particles
                potentialTemp.segment(count, tree[j][k]->potential.size()) = tree[j][k]->potential;
                count += tree[j][k]->potential.size();
            }
            potential += potentialTemp;
        }
    }

    Vec getTruePoten(Vec &charge)
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
            true_poten(i) = K->getRow(i, v).transpose() * charge;
        }
        return true_poten;
    }
    // Total relative error
    double compute_total_error()
    {
        Vec true_poten(N);
        std::vector<int> ind;
        for (int i = 0; i < N; ++i)
        {
            ind.push_back(i);
        }
        #pragma omp parallel for
        for (int i = 0; i < N; ++i)
        {
            true_poten(i) = K->getRow(i, ind).transpose() * collective_charge;
        }
        Vec errVec = true_poten - collective_potential;
        double error = errVec.norm() / true_poten.norm();
        return error;
    }
};

#endif