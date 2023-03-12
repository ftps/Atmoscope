#ifndef TEST_HPP
#define TEST_HPP

#include "turb.hpp"
#include <functional>
#include <algorithm>

void plotLegPoly(const int& l, const int& m);
void plotLs(const int& l, const int& ngraph = 100);
void plotLsSphere(const int& l, const int& ngraph = 100);
void plotSphereHarm(const int& l, const int& m, const int& N = 100);
void plotEnergy(const int& l, const int& n = 1);
void plotMap(const Turbulence& t);
void plotCn2();
void monteCarloCn2(const Turbulence& t, const uint& n = 100000);



std::array<double,2> linReg(const VpairXY& xy);
void testLinReg(const double& m = 1, const double& b = 0, const double& e = 0.1);

#endif