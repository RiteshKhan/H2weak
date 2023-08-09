#ifndef Integral2D_HPP
#define Integral2D_HPP
#include "Gauss_Legendre_Nodes_and_Weights.hpp"
#include "headers.hpp"

int n_gauss_nodes = 16;

#ifdef USE_real
dtype Integrand(double x, double y)
{
    double c = -1.0 / (2 * PI);
    return c * 0.5 * log(x * x + y * y);
}
#endif

#ifdef USE_Hankel

double besselJ(int n, double x)
{
    if (n >= 0)
    {
        double temp = boost::math::cyl_bessel_j(double(n), x);
        return temp;
    }
    else
    {
        double temp = boost::math::cyl_bessel_j(double(-n), x);
        if (-n % 2 == 0)
            return temp;
        else
            return -temp;
    }
}

double besselY(int n, double x)
{
    if (n >= 0)
    {
        double temp = boost::math::cyl_neumann(double(n), x);
        return temp;
    }
    else
    {
        double temp = boost::math::cyl_neumann(double(-n), x);
        if (-n % 2 == 0)
            return temp;
        else
            return -temp;
    }
}

dtype Integrand(double x, double y)
{
    double R = std::sqrt(x * x + y * y);
    return I * (besselJ(0, R) + I * besselY(0, R)) / 4.0;
}
#endif

// a = [x_low,y_low] and b = [x_upper,y_upper]
dtype double_integral(double *&a, double *&b)
{
    dtype result = 0.0;
    double tx, ty, cx, cy, Lx, Ly;
    double *nodes_x, *weights_x;
    Gauss_Legendre_Nodes_and_Weights(n_gauss_nodes, nodes_x, weights_x);
    double *nodes_y, *weights_y;
    Gauss_Legendre_Nodes_and_Weights(n_gauss_nodes, nodes_y, weights_y);
    // Shift the nodes from [-1,1] to that interval
    cx = 0.5 * (b[0] - a[0]);
    cy = 0.5 * (b[1] - a[1]);

    Lx = 0.5 * (b[0] + a[0]); // half length
    Ly = 0.5 * (b[1] + a[1]);

    // Gauss-Legendre Quadrature formula
    for (int i = 0; i < n_gauss_nodes; i++)
    {
        tx = cx * nodes_x[i] + Lx;
        for (int j = 0; j < n_gauss_nodes; j++)
        {
            ty = cy * nodes_y[j] + Ly;
            result += weights_x[i] * weights_y[j] * Integrand(tx, ty);
        }
    }
    // Final result due to change of variables
    result *= (cx * cy);
    result *= 4; // This is reducing the integral into half domain
    // std::cout << "The double integral is " << result << std::endl;
    return result;
}
#endif