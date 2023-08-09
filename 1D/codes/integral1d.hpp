#ifndef Integral1D_HPP
#define Integral1D_HPP
#include "Gauss_Legendre_Nodes_and_Weights.hpp"
#include "headers.hpp"

int n_gauss_nodes = 16;

#ifdef USE_real
dtype Integrand(double x)
{
    return log(fabs(1.0 * x));
}
#endif


// a = [x_low,y_low] and b = [x_upper,y_upper]
dtype single_integral(double *&a, double *&b)
{
    dtype result = 0.0;
    double tx, cx, Lx;
    double *nodes_x, *weights_x;
    Gauss_Legendre_Nodes_and_Weights(n_gauss_nodes, nodes_x, weights_x);
    // Shift the nodes from [-1,1] to that interval
    cx = 0.5 * (b[0] - a[0]);

    Lx = 0.5 * (b[0] + a[0]); // half length

    // Gauss-Legendre Quadrature formula
    for (int i = 0; i < n_gauss_nodes; i++)
    {
        tx = cx * nodes_x[i] + Lx;
        result += weights_x[i] * Integrand(tx);
    }
    // Final result due to change of variables
    result *= cx;
    result *= 2; // This is reducing the integral into half domain
    // std::cout << "The double integral is " << result << std::endl;
    return result;
}
#endif