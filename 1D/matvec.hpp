#ifndef __INITIALISATION_HPP__
#define __INITIALISATION_HPP__

#include "HODLR1Dtree.hpp"

class initialisation
{
public:
    std::vector<double> particles;
    int nParticlesInLeaf;
    double L;
    int TOL_POW;
    int Kernel_Choice;
    int N;
    int nLevels;
    hodlr1d_tree<userkernel> *K = NULL;

    initialisation(int N, int nLevels, int nParticlesInLeaf, double L, int TOL_POW, int Kernel_Choice)
    {
        this->nParticlesInLeaf = nParticlesInLeaf;
        this->L = L;
        this->TOL_POW = TOL_POW;
        this->N = N; // Number of particles.
        this->nLevels = nLevels;
        this->Kernel_Choice = Kernel_Choice;
    }
    ~initialisation(){};

#ifdef USE_nHODLRdD
    void tree_initialisation()
    {
        userkernel *mykernel = new userkernel(particles,Kernel_Choice);
        this->K = new hodlr1d_tree<userkernel>(mykernel, N, nLevels, nParticlesInLeaf, L, TOL_POW);

        K->set_Uniform_Nodes();

        this->particles = K->K->particles;
        this->K->K = new userkernel(particles, Kernel_Choice);
        K->set_Uniform_Nodes();
        K->createTree();
        K->assign_all_interaction();
        K->assign_center_location();

        K->assignChargeLocations();
        K->assignNonLeafChargeLocations();

        K->getNodes();

        K->assemble_M2M();
        K->assemble_M2L();
        K->assemble_NearField();
    }

    void mat_vec_prod(Vec &charge, Vec &output)
    {
        K->assignCharges(charge);

        K->evaluate_M2M();
        K->evaluate_M2L();
        K->evaluate_L2L();

        K->evaluate_NearField();

        K->collect_all_potentials(output);
        K->reorder(output);
        // std::cout << output << std::endl;
        std::cout << "Assembly time for MAT-VEC: " << K->elapsed_assem << std::endl;
        std::cout << "MAT-VEC time: " << K->elapsed_mvp << std::endl;
        std::cout << "The error is: " << K->compute_total_error() << std::endl;
        std::cout << std::endl;
        std::cout << "===========================" << std::endl;
    }
#endif

#ifdef USE_snHODLRdD
    void tree_initialisation()
    {
        userkernel *mykernel = new userkernel(particles,Kernel_Choice);
        this->K = new hodlr1d_tree<userkernel>(mykernel, N, nLevels, nParticlesInLeaf, L, TOL_POW);

        K->set_Uniform_Nodes();

        this->particles = K->K->particles;
        this->K->K = new userkernel(particles, Kernel_Choice);
        K->set_Uniform_Nodes();
        K->createTree();
        K->assign_all_interaction();
        K->assign_center_location();

        K->assignChargeLocations();
        K->assignNonLeafChargeLocations();

        K->getUVtTree_nn();
        K->assemble_NearField();
    }

    void mat_vec_prod(Vec &charge, Vec &output)
    {
        K->assignCharges(charge);

        K->fast_multiplication_nn();
        K->evaluate_NearField();
        

        K->collect_all_potential_nn(output);
        K->reorder(output);
        std::cout << "Assembly time for MAT-VEC: " << K->elapsed_assem << std::endl;
        std::cout << "MAT-VEC time: " << K->elapsed_mvp << std::endl;
        std::cout << "The error is: " << K->compute_total_error() << std::endl;
        std::cout << std::endl;
        std::cout << "===========================" << std::endl;
    }
#endif

};
#endif