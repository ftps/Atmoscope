#include "poly.hpp"

long double factTR(const long double& n, const long double& res)
{
    if(n <= 1) return res;
    return factTR(n-1, res*n);
}

long double fact(const int& n)
{
    return factTR(n-1, (n) ? n : 1);
}

long double ring(const int& l)
{
    long double aux = 1;
    
    for(int k = 0; k < l; ++k){
        aux *= (l - k - 0.5 );
    }

    return aux;
}




Poly::Poly(const int& l, const int& m) : l(l), m(m) { }

long double Poly::operator()(const double& x) const
{
    uint i = std::floor(((x + 1.0)/2.0)*(vals.size()-1));
    if(i == vals.size()-1) return vals.back();

    return std::lerp(vals[i], vals[i+1], (x+1 - i*dx)/dx);
    //return (vals[i+1] - vals[i])*(x+1 - i*dx)/dx + vals[i];
}

void Poly::addVal(const long double& y)
{
    vals.emplace_back(y);
    dx = 2.0/(vals.size()-1);
}

VpairXY Poly::getVals() const
{
    VpairXY res;

    for(uint i = 0; i < vals.size(); ++i){
        res.emplace_back(i*dx-1, vals[i]);
    }

    return res;
}

VpairXYld Poly::getValsld() const
{
    VpairXYld res;

    for(uint i = 0; i < vals.size(); ++i){
        res.emplace_back(i*dx-1, vals[i]);
    }

    return res;
}





long double PolyCalc::calcPoly(const int& l, const int& m, const double& x)
{
    long double aux;

    if(m > l) return 0;
    else if(l == m){
        aux = powq(-2.0, l)*powq(std::abs(1.0 - x*x), l/2.0)*ring(l);
    }
    else{
        aux = calcPoly(l-1, m, x); // = res[m][l-1]
        aux = ((2*l-1)*x*aux - ((l > 1) ? (l+m-1)*res[m][l-2] : 0.0))/(l-m);
    }

    res[m][l] = aux;
    return aux;
}


void PolyCalc::calcAll(const int& lmax, const int& N)
{
    double dx = 2.0/(N-1);
    int n = 0;

    this->lmax = lmax;
    this->N = N;

    res.clear();
    res.resize(lmax+1);
    poly.clear();
    poly.resize(lmax+1);
    normCoef.clear();
    //normCoef.resize(lmax+1);

    for(std::vector<long double>& vlm : res){
        vlm.clear();
        vlm.resize(lmax+1, 0);
    }
    for(int l = 0; l <= lmax; ++l){
        normCoef.emplace_back(std::vector<double>(l+1));
        for(int m = 0; m <= l; ++m){
            poly[l].push_back(Poly(l, m));
            if(l) normCoef.back()[m] = (m) ? sqrtq(2*((2*l+1)/(4*M_PI)) * fact(l - m)/fact(l+m)) : sqrtq((2*l+1)/(4*M_PI));
        }
    }

    for(double x = -1; x < 1+0.1*dx; x += dx){
        if((n++)%50 == 0) std::cout << "At iteration " << n << std::endl;
        for(int m = 0; m <= lmax; ++m){
            calcPoly(lmax, m, x);
            for(int l = m; l <= lmax; ++l){
                poly[l][m].addVal(res[m][l]);
            }
        }
    }
}

void PolyCalc::printPoly()
{
    std::fstream fp;
    char buff[128];

    for(int l = 1; l <= lmax; ++l){
        fp.open("poly/P" + std::to_string(l) + ".dat", std::ofstream::out);
        fp << N << std::endl;
        for(const Poly& p : poly[l]){
            for(const pairXYld& xy : p.getValsld()){
                quadmath_snprintf(buff, 128, "%.6Qe", xy.second);
                fp << std::string(buff, std::find(buff, buff+128, '\0')) << " ";
            }
            fp << std::endl;
        }
        fp.close();
    }
}