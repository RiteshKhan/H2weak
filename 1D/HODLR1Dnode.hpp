#ifndef NODE_HPP
#define NODE_HPP

#include "headers.hpp"
class hodlr_node{
    public:
        int nodeNumber;
        int parentNumber;
        int* childrenNumbers = new int[2];
        int* neighborNumbers = new int[2];
        double center;
        Vec charges, potential, potential_nn;
        std::vector<int> chargeLocations;
        std::vector<int> incoming_checkPoints, incoming_chargePoints;
        Vec outgoing_charges;
        Vec incoming_potential;
        std::map<int, Mat> M2L;
        Mat denseMatrices;
        Mat *M2M = new Mat[2];
        Mat L, R;
        std::map<int, Mat> U, V;
        hodlr_node(){
            this->nodeNumber   = -1;
            this->parentNumber = -1;
            for(int i=0; i<2; ++i){
                childrenNumbers[i] = -1;
                neighborNumbers[i] = -1;
            }
        }
        ~hodlr_node() {};
};

#endif