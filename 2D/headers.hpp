#ifndef Headers_HPP
#define Headers_HPP

#include <iostream>
#include <set>
#include <vector>
#include <Eigen/Dense>
#include <map>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <ios>

const double PI = 4.0 * atan(1);
const std::complex<double> I(0.0, 1.0);

struct pts2D
{
	double x, y;
};

double dist2D(pts2D &p1, pts2D &p2)
{
    double X = (double)(p1.x - p2.x) * (p1.x - p2.x);
    double Y = (double)(p1.y - p2.y) * (p1.y - p2.y);
    return (double)sqrt(X + Y);
}


#ifdef USE_Hankel
using dtype = std::complex<double>;
using dtype_base = double;
using Mat = Eigen::MatrixXcd;
using Vec = Eigen::VectorXcd;
#include <boost/math/special_functions/bessel.hpp>
#endif

#ifdef USE_real
using dtype = double;
using dtype_base = double;
using Mat = Eigen::MatrixXd;
using Vec = Eigen::VectorXd;
#endif

#endif