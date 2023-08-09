//
// Remodified from : https://github.com/vaishna77/NNCA2D/FMM2DTree.hpp to keep comparable with FMM or H matrix
// HODLR2D cluster
// 
#ifndef __HODLR2DNODE_HPP__
#define __HODLR2DNODE_HPP__
#include "headers.hpp"
class HODLR2DNode
{
public:
    int boxNumber;
    int parentNumber;
    int *childrenNumbers = new int[4];
    int *neighborNumbers = new int[4];
    std::vector<int> interactionList; // To store the interaction list of a node
    std::map<int, bool> isVertex;       // Tag vertex or non-vertex. It will help to separate out the vertex-sharing and far-field interactions.

    std::vector<int> chargeLocations;
    std::vector<int> incoming_checkPoints, incoming_checkPoints_ver, incoming_chargePoints, incoming_chargePoints_ver;
    pts2D center;
    Vec charges;
    Vec outgoing_charges;
    Vec incoming_potential, incoming_potential_ver;
    Vec potential, potential_ver;
    Mat *L2P = new Mat[1];
    Mat *L_ver = new Mat[1];
    Mat *R_ver = new Mat[1];
    Mat *M2M_ver = new Mat[4]; // M2M for vertex sharing interaction
    Mat *denseMatrices = new Mat[5];
    std::map<int, Mat> M2L_far, M2L_ver;
    std::map<int, Mat> Ac, Ar;
    Mat *L = new Mat[36]; 
    Mat *R = new Mat[36];
    std::vector<int> *row_basis = new std::vector<int>[36];
    std::vector<int> *col_basis = new std::vector<int>[36];

    // Box constructor
    HODLR2DNode()
    {
        this->boxNumber = -1;
        this->parentNumber = -1;
        for (int l = 0; l < 4; ++l)
        {
            this->childrenNumbers[l] = -1;
        }
        for (int l = 0; l < 4; ++l)
        {
            this->neighborNumbers[l] = -1;
        }
    }
    ~HODLR2DNode(){};
};
#endif