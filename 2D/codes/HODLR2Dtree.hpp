//
// Remodified from : https://github.com/vaishna77/NNCA2D/FMM2DTree.hpp to keep comparable with FMM or H matrix
// HODLR2D Tree
// 
#ifndef __HODLR2DTREE_HPP__
#define __HODLR2DTREE_HPP__
#include "aca.hpp"
#include "kernel.hpp"
#include "HODLR2Dnode.hpp"
// Tree class
template <typename kerneltype>
class HODLR2DTree
{
public:
    kerneltype *K;
    int nLevels;            //	Number of levels in the tree.
    int N;                  //	Number of particles.
    double L;               //	Semi-length of the simulation box.

    std::vector<int> nBoxesPerLevel;           //	Number of boxes at each level in the tree.
    std::vector<double> boxRadius;             //	Box radius at each level in the tree assuming the box at the root is [-1,1]^2
    std::vector<std::vector<HODLR2DNode*> > tree; //	The tree storing all the information.

    int nParticlesInLeafAlong1D;
    int nParticlesInLeaf;
    std::vector<double> Nodes1D;
    std::vector<pts2D> Nodes;
    std::vector<pts2D> gridPoints; // all particles in domain
    int TOL_POW;
    int sqrtRootN;
    Vec collective_charge;
    Vec collective_potential;

    double elapsed_mvp = 0.0;
    double elapsed_assem = 0.0;
    double memory = 0.0;
    HODLR2DTree(kerneltype *K, int sqrtRootN, int nLevels, int nParticlesInLeafAlong1D, double L, int TOL_POW)
    {
        this->K = K;
        this->sqrtRootN = sqrtRootN;
        this->nLevels = nLevels;
        this->L = L;
        this->nParticlesInLeafAlong1D = nParticlesInLeafAlong1D;
        this->nParticlesInLeaf = nParticlesInLeafAlong1D * nParticlesInLeafAlong1D;
        this->TOL_POW = TOL_POW;
        nBoxesPerLevel.push_back(1);
        boxRadius.push_back(L);
        for (int k = 1; k <= nLevels; ++k)
        {
            nBoxesPerLevel.push_back(4 * nBoxesPerLevel[k - 1]);
            boxRadius.push_back(0.5 * boxRadius[k - 1]);
        }
        this->N = sqrtRootN * sqrtRootN;
    }
    ~HODLR2DTree(){};
    // Set nodes. It will be location
    void set_Standard_Cheb_Nodes()
    {
        for (int k = 0; k < sqrtRootN; ++k)
        {
            Nodes1D.push_back(-cos((k + 0.5) / sqrtRootN * PI));
        }
        pts2D temp1;
        for (int j = 0; j < sqrtRootN; ++j)
        {
            for (int k = 0; k < sqrtRootN; ++k)
            {
                temp1.x = Nodes1D[k];
                temp1.y = Nodes1D[j];
                K->particles.push_back(temp1);
            }
        }
    }

    void set_Uniform_Nodes()
    {
        for (int k = 0; k < sqrtRootN; ++k)
        {
            Nodes1D.push_back(-L + 2.0 * L * (k + 1.0) / (sqrtRootN + 1.0));
        }
        pts2D temp1;
        for (int j = 0; j < sqrtRootN; ++j)
        {
            for (int k = 0; k < sqrtRootN; ++k)
            {
                temp1.x = Nodes1D[k];
                temp1.y = Nodes1D[j];
                K->particles.push_back(temp1);
            }
        }
    }
    // Need to shift location at leaf level
    void shift_Nodes(double radius, pts2D center, std::vector<pts2D> &particle_loc)
    {
        for (size_t i = 0; i < Nodes.size(); i++)
        {
            pts2D temp;
            temp.x = Nodes[i].x * radius + center.x;
            temp.y = Nodes[i].y * radius + center.y;
            particle_loc.push_back(temp);
        }
    }
    // Create tree and interaction list
    void createTree()
    {
        //	First create root and add to tree
        HODLR2DNode* root = new HODLR2DNode;
        root->boxNumber = 0;
        root->parentNumber = -1;
        #pragma omp parallel for
        for (int l = 0; l < 4; ++l)
        {
            root->childrenNumbers[l] = l;
        }
        #pragma omp parallel for
        for (int l = 0; l < 4; ++l)
        {
            root->neighborNumbers[l] = -1;
        }
        std::vector<HODLR2DNode*> rootLevel;
        rootLevel.push_back(root);
        tree.push_back(rootLevel);

        for (int j = 1; j <= nLevels; ++j)
        {
            std::vector<HODLR2DNode*> level;
            for (int k = 0; k < nBoxesPerLevel[j]; ++k)
            {
                HODLR2DNode* temp = new HODLR2DNode;
                temp->boxNumber = k;
                temp->parentNumber = k / 4;
                for (int l = 0; l < 4; ++l)
                {
                    temp->childrenNumbers[l] = 4 * k + l;
                }
                level.push_back(temp);
            }
            tree.push_back(level);
        }
    }
    //	Assigns the interactions for child0 of a box
    void assign_Child0_Interaction(int j, int k)
    {
        int nL = j + 1;
        int nC = 4 * k;
        int nN;

        //	Assign siblings
        {
            tree[nL][nC]->interactionList.push_back(nC + 2);
            tree[nL][nC]->neighborNumbers[1] = nC + 1;
            tree[nL][nC]->neighborNumbers[2] = nC + 3;
        }

        //	Assign children of parent's zeroth neighbor
        {

            nN = tree[j][k]->neighborNumbers[0];
            if (nN != -1)
            {
                tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[0]);
                tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[1]);
                tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[2]);
                tree[nL][nC]->neighborNumbers[0] = tree[j][nN]->childrenNumbers[3];
            }
        }

        //	Assign children of parent's first neighbor
        {
            nN = tree[j][k]->neighborNumbers[1];
            if (nN != -1)
            {
                tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[0]);
                tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[1]);
                tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[2]);
                tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[3]);
            }
        }

        //	Assign children of parent's second neighbor
        {

            nN = tree[j][k]->neighborNumbers[2];
            if (nN != -1)
            {
                tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[0]);
                tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[1]);
                tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[2]);
                tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[3]);
            }
        }

        //	Assign children of parent's third neighbor
        {

            nN = tree[j][k]->neighborNumbers[3];
            if (nN != -1)
            {
                tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[0]);
                tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[2]);
                tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[3]);
                tree[nL][nC]->neighborNumbers[3] = tree[j][nN]->childrenNumbers[1];
            }
        }
    }

    //	Assigns the interactions for child1 of a box
    void assign_Child1_Interaction(int j, int k)
    {
        int nL = j + 1;
        int nC = 4 * k + 1;
        int nN;

        //	Assign siblings
        {
            tree[nL][nC]->interactionList.push_back(nC + 2);
            tree[nL][nC]->neighborNumbers[3] = nC - 1;
            tree[nL][nC]->neighborNumbers[2] = nC + 1;
        }

        //	Assign children of parent's zeroth neighbor
        {

            nN = tree[j][k]->neighborNumbers[0];
            if (nN != -1)
            {
                tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[0]);
                tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[1]);
                tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[3]);
                tree[nL][nC]->neighborNumbers[0] = tree[j][nN]->childrenNumbers[2];
            }
        }

        //	Assign children of parent's first neighbor
        {

            nN = tree[j][k]->neighborNumbers[1];
            if (nN != -1)
            {
                tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[1]);
                tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[2]);
                tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[3]);
                tree[nL][nC]->neighborNumbers[1] = tree[j][nN]->childrenNumbers[0];
            }
        }

        //	Assign children of parent's second neighbor
        {

            nN = tree[j][k]->neighborNumbers[2];
            if (nN != -1)
            {
                tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[0]);
                tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[1]);
                tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[2]);
                tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[3]);
            }
        }

        //	Assign children of parent's third neighbor
        {

            nN = tree[j][k]->neighborNumbers[3];
            if (nN != -1)
            {
                tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[0]);
                tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[1]);
                tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[2]);
                tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[3]);
            }
        }
    }

    //	Assigns the interactions for child2 of a box
    void assign_Child2_Interaction(int j, int k)
    {
        int nL = j + 1;
        int nC = 4 * k + 2;
        int nN;

        //	Assign siblings
        {

            tree[nL][nC]->interactionList.push_back(nC - 2);
            tree[nL][nC]->neighborNumbers[0] = nC - 1;
            tree[nL][nC]->neighborNumbers[3] = nC + 1;
        }

        //	Assign children of parent's zeroth neighbor
        {

            nN = tree[j][k]->neighborNumbers[0];
            if (nN != -1)
            {
                tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[0]);
                tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[1]);
                tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[2]);
                tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[3]);
            }
        }

        //	Assign children of parent's first neighbor
        {

            nN = tree[j][k]->neighborNumbers[1];
            if (nN != -1)
            {
                tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[0]);
                tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[1]);
                tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[2]);
                tree[nL][nC]->neighborNumbers[1] = tree[j][nN]->childrenNumbers[3];
            }
        }

        //	Assign children of parent's second neighbor
        {

            nN = tree[j][k]->neighborNumbers[2];
            if (nN != -1)
            {
                tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[0]);
                tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[2]);
                tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[3]);
                tree[nL][nC]->neighborNumbers[2] = tree[j][nN]->childrenNumbers[1];
            }
        }

        //	Assign children of parent's third neighbor
        {

            nN = tree[j][k]->neighborNumbers[3];
            if (nN != -1)
            {
                tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[0]);
                tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[1]);
                tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[2]);
                tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[3]);
            }
        }
    }

    //	Assigns the interactions for child3 of a box
    void assign_Child3_Interaction(int j, int k)
    {
        int nL = j + 1;
        int nC = 4 * k + 3;
        int nN;

        //	Assign siblings
        {
            tree[nL][nC]->interactionList.push_back(nC - 2);
            tree[nL][nC]->neighborNumbers[0] = nC - 3;
            tree[nL][nC]->neighborNumbers[1] = nC - 1;
        }

        //	Assign children of parent's zeroth neighbor
        {

            nN = tree[j][k]->neighborNumbers[0];
            if (nN != -1)
            {
                tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[0]);
                tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[1]);
                tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[2]);
                tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[3]);
            }
        }

        //	Assign children of parent's first neighbor
        {

            nN = tree[j][k]->neighborNumbers[1];
            if (nN != -1)
            {
                tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[0]);
                tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[1]);
                tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[2]);
                tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[3]);
            }
        }

        //	Assign children of parent's second neighbor
        {

            nN = tree[j][k]->neighborNumbers[2];
            if (nN != -1)
            {
                tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[1]);
                tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[2]);
                tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[3]);
                tree[nL][nC]->neighborNumbers[2] = tree[j][nN]->childrenNumbers[0];
            }
        }

        //	Assign children of parent's third neighbor
        {

            nN = tree[j][k]->neighborNumbers[3];
            if (nN != -1)
            {
                tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[0]);
                tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[1]);
                tree[nL][nC]->interactionList.push_back(tree[j][nN]->childrenNumbers[3]);
                tree[nL][nC]->neighborNumbers[3] = tree[j][nN]->childrenNumbers[2];
            }
        }
    }

    // Interaction of all levels accross the tree
    void assign_all_interaction()
    {
        for (int j = 0; j < nLevels; ++j)
        {
            #pragma omp parallel for
            for (int k = 0; k < nBoxesPerLevel[j]; ++k)
            {
                assign_Child0_Interaction(j, k);
                assign_Child1_Interaction(j, k);
                assign_Child2_Interaction(j, k);
                assign_Child3_Interaction(j, k);
            }
        }
    }

    // Center location of different boxes
    void assign_center_location()
    {
        int J;
        tree[0][0]->center.x = 0.0;
        tree[0][0]->center.y = 0.0;
        for (int j = 0; j < nLevels; ++j)
        {
            J = j + 1;
            double shift = 0.5 * boxRadius[j];
            #pragma omp parallel for
            for (int k = 0; k < nBoxesPerLevel[j]; ++k)
            {
                tree[J][4 * k]->center.x = tree[j][k]->center.x - shift;
                tree[J][4 * k + 1]->center.x = tree[j][k]->center.x + shift;
                tree[J][4 * k + 2]->center.x = tree[j][k]->center.x + shift;
                tree[J][4 * k + 3]->center.x = tree[j][k]->center.x - shift;

                tree[J][4 * k]->center.y = tree[j][k]->center.y - shift;
                tree[J][4 * k + 1]->center.y = tree[j][k]->center.y - shift;
                tree[J][4 * k + 2]->center.y = tree[j][k]->center.y + shift;
                tree[J][4 * k + 3]->center.y = tree[j][k]->center.y + shift;
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

                    double temp = dist2D(tree[j][k]->center, tree[j][ki]->center);

                    if (temp - (2.0 * sqrt(2) * boxRadius[j]) <= 1.0e-16)
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

    // Assign charge location
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
                int Kp = 4 * k;
                for (size_t i = 0; i < tree[j][k]->chargeLocations.size(); i++)
                {
                    int index = tree[j][k]->chargeLocations[i];
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
            }
        }
    }

    // Assign charge at non leaf
    void assignNonLeafChargeLocations()
    {
        for (int j = nLevels - 1; j > 0; --j)
        {
            for (int k = 0; k < nBoxesPerLevel[j]; ++k)
            {
                tree[j][k]->chargeLocations.clear();
                for (int c = 0; c < 4; ++c)
                {
                    tree[j][k]->chargeLocations.insert(tree[j][k]->chargeLocations.end(), tree[j + 1][4 * k + c]->chargeLocations.begin(), tree[j + 1][4 * k + c]->chargeLocations.end());
                }
            }
        }
    }

    // Bottom to top pivot selection. We choose charge location as pivot only at leaf level.
    void pivot_selection_bottom_up(int j, int k, std::vector<int> &searchNodes)
    {
        if (j == nLevels)
        {
            searchNodes.insert(searchNodes.end(), tree[j][k]->chargeLocations.begin(), tree[j][k]->chargeLocations.end());
        }
        else
        {
            for (int c = 0; c < 4; ++c)
            {
                searchNodes.insert(searchNodes.end(), tree[j + 1][4 * k + c]->incoming_checkPoints.begin(), tree[j + 1][4 * k + c]->incoming_checkPoints.end());
            }
        }
    }

	void getNodes_incoming_box(int j, int k, int &ComputedRank)
	{
		std::vector<int> boxA_Nodes;
		pivot_selection_bottom_up(j, k, boxA_Nodes);
		std::vector<int> IL_Nodes; // indices
		for (size_t in = 0; in < tree[j][k]->interactionList.size(); ++in)
		{
			int kIL = tree[j][k]->interactionList[in];
			if (!tree[j][k]->isVertex[kIL])
			{
				pivot_selection_bottom_up(j, kIL, IL_Nodes);
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

    // Get nodes accross all levels
    void getNodes()
    {
        for (int j = nLevels; j >= 2; --j)
        {
            int ComputedRank;
            for (int k = 0; k < nBoxesPerLevel[j]; ++k)
            {
                getNodes_incoming_box(j, k, ComputedRank);
            }
        }
    }

    // Assign charges
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

    // Evaluate M2M bottom to top traverse
    void evaluate_M2M_far()
    {
        for (int j = nLevels; j >= 2; --j)
        {
            #pragma omp parallel for
            for (int k = 0; k < nBoxesPerLevel[j]; ++k)
            {
                Vec source_density;
                if (j == nLevels)
                {
                    source_density = tree[j][k]->charges;
                }
                else
                {
                    int vec_len = 0;
                    for (int c = 0; c < 4; ++c)
                    {
                        vec_len += tree[j + 1][4 * k + c]->outgoing_charges.size();
                    }
                    source_density = Vec::Zero(vec_len);
                    int ind = 0;
                    for (int c = 0; c < 4; ++c)
                    {
                        int l = tree[j + 1][4 * k + c]->outgoing_charges.size();
                        source_density.segment(ind, l) = tree[j + 1][4 * k + c]->outgoing_charges;
                        ind += l;
                    }
                }

                clock_t start;
                start = clock();
                tree[j][k]->outgoing_charges = tree[j][k]->L2P[0].transpose() * source_density;
                double duration = (clock() - start) / (double)CLOCKS_PER_SEC;
                elapsed_mvp += duration;
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

    // Evaluate M2L accross all levels
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

    // Evaluate L2L top to bottom traverse
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
                    for (int c = 0; c < 4; ++c)
                    {
                        tree[j + 1][4 * k + c]->incoming_potential += temp.segment(ind, tree[j + 1][4 * k + c]->incoming_checkPoints.size());
                        ind += tree[j + 1][4 * k + c]->incoming_checkPoints.size();
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

    ///////////////////////////////////////////////////////////////////
    /*********** potential for vertex sharng interaction *************/
    ///////////////////////////////////////////////////////////////////

    // Top-bottom approach
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
			IL_Nodes.insert(IL_Nodes.end(), tree[j - 1][k / 4]->incoming_chargePoints_ver.begin(), tree[j - 1][k / 4]->incoming_chargePoints_ver.end());
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
        for (int j = 1; j <= nLevels; j++)
        {
            for (int k = 0; k < nBoxesPerLevel[j]; k++)
            {
                if (j == nLevels)
                {
                    Mat temp = K->getMatrix(tree[j][k]->incoming_chargePoints_ver, tree[j][k]->chargeLocations);
                    Mat temp2 = tree[j][k]->R_ver[0].transpose().triangularView<Eigen::Lower>().solve(temp);
                    tree[j][k]->M2M_ver[0] = tree[j][k]->L_ver[0].transpose().triangularView<Eigen::Upper>().solve(temp2);
                }
                else
                {
                    for (int c = 0; c < 4; c++)
                    {
                        Mat temp = K->getMatrix(tree[j][k]->incoming_chargePoints_ver, tree[j + 1][4 * k + c]->incoming_checkPoints_ver);
                        Mat temp2 = tree[j][k]->R_ver[0].transpose().triangularView<Eigen::Lower>().solve(temp);
                        tree[j][k]->M2M_ver[c] = tree[j][k]->L_ver[0].transpose().triangularView<Eigen::Upper>().solve(temp2);
                    }
                }
            }
        }
    }

    void evaluate_M2M_ver()
    {
        for (int j = nLevels; j >= 1; j--)
        {
            #pragma omp parallel for
            for (int k = 0; k < nBoxesPerLevel[j]; k++)
            {
                if (j == nLevels)
                {
                    clock_t start8;
                    start8 = clock();
                    Mat M2M = tree[j][k]->M2M_ver[0];
                    Vec w = tree[j][k]->charges;
                    elapsed_assem += (clock() - start8) / (double)CLOCKS_PER_SEC;


                    clock_t start;
                    start = clock();
                    tree[j][k]->outgoing_charges = M2M * w;
                    elapsed_mvp += (clock() - start) / (double)CLOCKS_PER_SEC;
                }
                else
                {
                    tree[j][k]->outgoing_charges = Vec::Zero(tree[j][k]->incoming_checkPoints_ver.size());
                    for (int c = 0; c < 4; c++)
                    {
                        clock_t start9;
                        start9 = clock();
                        Mat M2M = tree[j][k]->M2M_ver[c];
                        Vec w = tree[j + 1][4 * k + c]->outgoing_charges;
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
        for (int j = 1; j <= nLevels; j++)
        {
            #pragma omp parallel for
            for (int k = 0; k < nBoxesPerLevel[j]; k++)
            {
                if (j != nLevels)
                {
                    for (int c = 0; c < 4; c++)
                    {
                        clock_t start;
                        start = clock();
                        tree[j + 1][4 * k + c]->incoming_potential_ver += tree[j][k]->M2M_ver[c].transpose() * tree[j][k]->incoming_potential_ver;
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

    ////////////// Potential corresponding to vertices (ACA) //////////////

    void getUVtInstance(int j, int k, int ki, Mat &L, Mat &R, std::vector<int> &row_indices, std::vector<int> &col_indices)
    {
        int computed_rank;
        std::vector<int> row_indices_local, col_indices_local;
        LowRank *LR = new LowRank(K, TOL_POW, tree[j][k]->chargeLocations, tree[j][ki]->chargeLocations);
        LR->ACA_only_nodesCUR(row_indices_local, col_indices_local, computed_rank, L, R);
        for (int r = 0; r < computed_rank; r++)
        {
            row_indices.push_back(tree[j][k]->chargeLocations[row_indices_local[r]]);
        }
        for (int c = 0; c < computed_rank; c++)
        {
            col_indices.push_back(tree[j][ki]->chargeLocations[col_indices_local[c]]);
        }
        delete LR;
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
			for (int n = 0; n < 4; n++)
			{
				int nn = tree[nLevels][k]->neighborNumbers[n];
				if (nn != -1)
				{
					tree[nLevels][k]->denseMatrices[n] = K->getMatrix(tree[nLevels][k]->chargeLocations, tree[nLevels][nn]->chargeLocations);;
				}
			}
			tree[nLevels][k]->denseMatrices[4] = K->getMatrix(tree[nLevels][k]->chargeLocations, tree[nLevels][k]->chargeLocations);;
		}
	}

    // Evaluate leaf level computations, neighbor & self interaction
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
        for (int k = 0; k < nBoxesPerLevel[nLevels]; ++k)
        {
            for (int n = 0; n < 4; ++n)
            {
                int nn = tree[nLevels][k]->neighborNumbers[n];
                if (nn != -1)
                {
                    clock_t start;
                    start = clock();
                    tree[nLevels][k]->potential += tree[nLevels][k]->denseMatrices[n] * tree[nLevels][nn]->charges;
                    elapsed_mvp += (clock() - start) / (double)CLOCKS_PER_SEC;
                }
            }
            // Self interaction
            clock_t start;
            start = clock();
            tree[nLevels][k]->potential += tree[nLevels][k]->denseMatrices[4] * tree[nLevels][k]->charges;
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
            {
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
        for (int k = 0; k < nBoxesPerLevel[nLevels]; ++k)
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

	void findMemory_n() {
        // M2M + L2L memory
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
					for (int c = 0; c < 4; c++)
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

        // Dense matrices
		for (int k = 0; k < nBoxesPerLevel[nLevels]; k++) {
			memory += tree[nLevels][k]->chargeLocations.size()*tree[nLevels][k]->chargeLocations.size(); //self
			for (size_t n = 0; n < 4; n++) {
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
			for (size_t n = 0; n < 4; n++) {
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
			for (size_t n = 0; n < 4; n++) {
				int nn = tree[nLevels][k]->neighborNumbers[n];
				if(nn != -1) {
					memory += tree[nLevels][k]->chargeLocations.size()*tree[nLevels][nn]->chargeLocations.size();
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


    // Total norm 2 error
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