#ifndef __INITIALISATION_HPP__
#define __INITIALISATION_HPP__

#include "HODLR3Dtree.hpp"

class initialisation
{
public:
    std::vector<pts3D> particles;
    int nLevels;
    int nParticlesInLeafAlong1D;
    int nParticlesInLeaf;
    double L;
    int TOL_POW;
    int Kernel_Choice;
    int N;
    int cubeRootN;
    HODLR3DTree<userkernel> *K = NULL;

    initialisation(int cubeRootN, int nLevels, int nParticlesInLeafAlong1D, int L, int TOL_POW, int Kernel_Choice)
    {
        this->cubeRootN = cubeRootN;
        this->nLevels = nLevels;
        this->nParticlesInLeafAlong1D = nParticlesInLeafAlong1D;
        this->nParticlesInLeaf = nParticlesInLeafAlong1D * nParticlesInLeafAlong1D * nParticlesInLeafAlong1D;
        this->L = L;
        this->TOL_POW = TOL_POW;
        this->Kernel_Choice = Kernel_Choice;
    }
    ~initialisation(){};

#ifdef USE_nHODLRdD
    void tree_initialisation()
    {
        userkernel* mykernel = new userkernel(particles, Kernel_Choice); 
        this->K = new HODLR3DTree<userkernel>(mykernel, cubeRootN, nLevels, nParticlesInLeafAlong1D, L, TOL_POW);

        K->set_Uniform_Nodes();
        
        this->particles = K->K->particles;
        this->K->K = new userkernel(particles, Kernel_Choice);


        K->createTree();
        K->assign_Tree_Interactions();
        K->assign_Center_Location();

        K->separate_vertex_interaction();

        K->assignChargeLocations();
        K->assignNonLeafChargeLocations();

        K->getNodes();
        K->getNodes_ver();

        K->assemble_M2L_far();
        K->assemble_M2L_ver();
        
        K->assemble_M2M_ver();
        K->assemble_NearField();
    }

    void mat_vec_prod(Vec &charge, Vec &output)
    {
        K->assignCharges(charge);
        K->evaluate_M2M_far();
        K->evaluate_M2L_far();
        K->evaluate_L2L_far();
        
        K->evaluate_M2M_ver();
        K->evaluate_M2L_ver();
        K->evaluate_L2L_ver();

        K->evaluate_NearField();

        K->collect_all_potentials_n(output);
        K->reorder(output);
        // std:: cout << "Assembly time for MAT-VEC: " << K->elapsed_assem << std::endl;
        // std:: cout << "MAT-VEC time: " << K->elapsed_mvp << std::endl;
        // //std:: cout << "The error is: " << K->compute_total_error() << std::endl;
        // std::cout << std::endl;
        // std::cout << "===========================" << std::endl;
    }
#endif
#ifdef USE_snHODLRdD
    void tree_initialisation()
    {
        userkernel* mykernel = new userkernel(particles, Kernel_Choice); 
        this->K = new HODLR3DTree<userkernel>(mykernel, cubeRootN, nLevels, nParticlesInLeafAlong1D, L, TOL_POW);

        K->set_Uniform_Nodes();
        
        this->particles = K->K->particles;
        this->K->K = new userkernel(particles, Kernel_Choice);

        K->createTree();
        K->assign_Tree_Interactions();
        K->assign_Center_Location();

        K->separate_vertex_interaction();


        K->assignChargeLocations();
        K->assignNonLeafChargeLocations();

        K->getNodes();
        K->assemble_M2L_far();
        K->assemble_NearField();

        K->getUVtTree_sn();
        K->assembleAcAr_sn();

    }

    void mat_vec_prod(Vec &charge, Vec &output)
    {
        K->assignCharges(charge);

        K->evaluate_M2M_far();
        K->evaluate_M2L_far();
        K->evaluate_L2L_far();

        K->evaluate_NearField();

        K->evaluateFarField_sn();


        K->collect_all_potentials_semi_1(output);
        K->collect_all_potentials_semi_2(output);
        K->reorder(output);
        // std:: cout << "Assembly time for MAT-VEC: " << K->elapsed_assem << std::endl;
        // std:: cout << "MAT-VEC time: " << K->elapsed_mvp << std::endl;
        // ///std:: cout << "The error is: " << K->compute_total_error() << std::endl;
        // std::cout << std::endl;
        // std::cout << "===========================" << std::endl;
    }
#endif
#ifdef USE_HODLRdD
    void tree_initialisation()
    {
        userkernel* mykernel = new userkernel(particles, Kernel_Choice); 
        this->K = new HODLR3DTree<userkernel>(mykernel, cubeRootN, nLevels, nParticlesInLeafAlong1D, L, TOL_POW);

        K->set_Uniform_Nodes();
        
        this->particles = K->K->particles;
        this->K->K = new userkernel(particles, Kernel_Choice);

        K->createTree();
        K->assign_Tree_Interactions();
        K->assign_Center_Location();


        K->assignChargeLocations();
        K->assignNonLeafChargeLocations();

        K->getUVtTree_nn();
        K->assembleAcAr_nn();

        K->assemble_NearField();
    }

    void mat_vec_prod(Vec &charge, Vec &output)
    {
        K->assignCharges(charge);

        K->evaluateFarField_nn();

        K->evaluate_NearField();

        K->collect_all_potentials_nn(output);
        K->reorder(output);
        // std:: cout << "Assembly time for MAT-VEC: " << K->elapsed_assem << std::endl;
        // std:: cout << "MAT-VEC time: " << K->elapsed_mvp << std::endl;
        // //std:: cout << "The error is: " << K->compute_total_error() << std::endl;
        // std::cout << std::endl;
        // std::cout << "===========================" << std::endl;
    }
#endif
};
#endif