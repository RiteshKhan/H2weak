//
// Remodified from : https://github.com/vaishna77/HODLR3D/ to keep comparable with FMM or H matrix
// HODLR3D cluster
//
#ifndef __HODLR3DNODE_HPP__
#define __HODLR3DNODE_HPP__
#include "headers.hpp"
class HODLR3DNode
{
public:
    int boxNumber;
    int parentNumber;
    int childrenNumbers[8];
    int neighborNumbers[26]; // other than self

    std::vector<int> interactionList;
    std::map<int, bool> isVertex;
    pts3D center;
    Vec charges, potential, potential_ver, collective_potential;
    Mat *L_far = new Mat[216]; // upper bound of neighbors for FMM3D is 27; upper bound of IL for FMM3D is 216=27*8;
    Mat *R_far = new Mat[216];

    std::vector<int> *row_basis_far = new std::vector<int>[216];
    std::vector<int> *col_basis_far = new std::vector<int>[216];
    std::vector<int> *row_basis_ver = new std::vector<int>[216];
    std::vector<int> *col_basis_ver = new std::vector<int>[216];

    Vec outgoing_charges;                                              // equivalent densities {f_{k}^{B,o}}
    Vec incoming_charges;                                              // equivalent densities {f_{k}^{B,i}}
    Vec incoming_potential, incoming_potential_ver;                    // check potentials {u_{k}^{B,i}}
    std::vector<int> incoming_chargePoints, incoming_chargePoints_ver; // equivalent points {y_{k}^{B,i}}
    std::vector<int> incoming_checkPoints, incoming_checkPoints_ver;   // check points {x_{k}^{B,i}}
    std::vector<int> chargeLocations;
    Mat *L2P = new Mat[1];
    Mat *L_ver = new Mat[1];
    Mat *R_ver = new Mat[1];
    Mat *M2M_ver = new Mat[8]; // M2M for vertex sharing interaction
    Mat *denseMatrices = new Mat[27];
    std::map<int, Mat> M2L_far, M2L_ver;
    std::map<int, Mat> Ac, Ar;
    std::vector<int> *row_basis = new std::vector<int>[216];
	std::vector<int> *col_basis = new std::vector<int>[216];
	Mat *L = new Mat[216]; // upper bound of neighbors for FMM3D is 27; upper bound of IL for FMM3D is 216=27*8;
	Mat *R = new Mat[216];
    HODLR3DNode()
    {
        this->boxNumber = -1;
        this->parentNumber = -1;
        for (int l = 0; l < 8; ++l)
        {
            this->childrenNumbers[l] = -1;
        }
        for (int l = 0; l < 26; ++l)
        {
            this->neighborNumbers[l] = -1;
        }
    }
    ~HODLR3DNode(){};
};
#endif