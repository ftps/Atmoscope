#ifndef SPHERE_HPP
#define SPHERE_HPP

#include "poly.hpp"
#include<filesystem>

namespace fs = std::filesystem;

class SphereHarm : public PolyCalc {
    public:
        SphereHarm();
        double operator()(const int& l, const int& m, const double& tet, const double& phi) const;
        double gphi(const int& l, const int& m, const double& tet, const double& phi) const;
        double gtet(const int& l, const int& m, const double& tet, const double& phi) const;

        //int lmax;
    //private:
        //std::vector<std::vector<Poly>> poly;
        //std::vector<std::vector<long double>> normCoef;
};

#endif