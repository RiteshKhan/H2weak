#ifndef Integral3D_HPP
#define Integral3D_HPP
#include "headers.hpp"
#include "Gauss_Legendre_Nodes_and_Weights.hpp"

int n_gauss_nodes = 16;
dtype Integrand(double x, double y, double z)
{
    return 1.0/((4*PI)*sqrt(x*x + y*y + z*z));
}

// a = [x_low,y_low] and b = [x_upper,y_upper]
dtype triple_integral(double*& a, double*& b){
    dtype result = 0.0;
    double tx, ty, tz, cx, cy, cz, Lx, Ly, Lz;
    double *nodes_x,*weights_x;
    Gauss_Legendre_Nodes_and_Weights(n_gauss_nodes, nodes_x, weights_x);
    double *nodes_y,*weights_y;
    Gauss_Legendre_Nodes_and_Weights(n_gauss_nodes, nodes_y, weights_y);
    double *nodes_z,*weights_z;
    Gauss_Legendre_Nodes_and_Weights(n_gauss_nodes, nodes_z, weights_z);
    // Shift the nodes from [-1,1] to that interval
    cx = 0.5 * (b[0] - a[0]);
    cy = 0.5 * (b[1] - a[1]);
    cz = 0.5 * (b[2] - a[2]);

    Lx = 0.5 * (b[0] + a[0]);   // half length
    Ly = 0.5 * (b[1] + a[1]);
    Lz = 0.5 * (b[2] + a[2]);

    // Gauss-Legendre Quadrature formula
    for(int i=0; i<n_gauss_nodes; i++)
    {
            tx = cx*nodes_x[i] + Lx;
        for(int j=0; j<n_gauss_nodes; j++)
        {
            ty = cy*nodes_y[j] + Ly;
            for(int k=0; k<n_gauss_nodes; k++)
            {
                tz = cz*nodes_z[k] + Lz;
                result += (weights_x[i]*weights_y[j]*weights_z[k] * Integrand(tx,ty,tz));
            }
        }
    }
    // Final result due to change of variables
    result *= (cx*cy*cz);
    result *= 8; // This is reducing the integral into half domain
    //std::cout << "The double integral is " << result << std::endl;
    return result;
}
#endif