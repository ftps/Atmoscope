#ifndef POLY_HPP
#define POLY_HPP

#include "gnuplot-iostream.hpp"
#include "utils.hpp"


long double fact(const int& n);
long double ring(const int& l);

class Poly {
    public:
        Poly(const int& l, const int& m);
        virtual long double operator()(const double& x) const;
        void addVal(const long double& y);
        VpairXY getVals() const;
        VpairXYld getValsld() const;

        const int l, m;
    private:
        double dx;
        std::vector<long double> vals;

};

class PolyCalc {
    public:
        void calcAll(const int& lmax = 250, const int& N = 100);
        void printPoly();
    //private:
        long double calcPoly(const int& l, const int& m, const double& x);

        int lmax = 0, N = 0;
        std::vector<std::vector<long double>> res;
        std::vector<std::vector<Poly>> poly;
        std::vector<std::vector<double>> normCoef;
};

#endif