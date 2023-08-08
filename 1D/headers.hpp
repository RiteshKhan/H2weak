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


#ifdef USE_Complex
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