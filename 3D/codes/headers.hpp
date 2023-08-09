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
struct pts3D
{
	double x, y, z;
};


// 3d distance, we need it to separate the vertex sharing interaction
double dist3D(pts3D &p1, pts3D &p2)
{
	double X = double((p1.x - p2.x) * (p1.x - p2.x));
	double Y = double((p1.y - p2.y) * (p1.y - p2.y));
	double Z = double((p1.z - p2.z) * (p1.z - p2.z));
	return sqrt(X + Y + Z);
}

#ifdef USE_Helm
using dtype = std::complex<double>;
using dtype_base = double;
using Mat = Eigen::MatrixXcd;
using Vec = Eigen::VectorXcd;
#endif

#ifdef USE_real
using dtype = double;
using dtype_base = double;
using Mat = Eigen::MatrixXd;
using Vec = Eigen::VectorXd;
#endif

#endif