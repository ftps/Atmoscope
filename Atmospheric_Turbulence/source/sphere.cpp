#include "sphere.hpp"

template<typename T>
T cotan(const T& x)
{
    return std::cos(x)/std::sin(x);
}

SphereHarm::SphereHarm()
{
    /*std::fstream fp;
    int N;
    double aux;
    std::string buf;
    poly.emplace_back(std::vector<Poly>());
    normCoef.emplace_back(std::vector<long double>());

    for(int l = 1; fs::exists(fs::path("poly/P" + std::to_string(l) + ".dat")); ++l){
        ++lmax;
        if(!(l%50)) std::cout << l << std::endl;
        poly.emplace_back(std::vector<Poly>());
        normCoef.emplace_back(std::vector<long double>());
        fp.open("poly/P" + std::to_string(l) + ".dat", std::ifstream::in);
        fp >> N;
        for(int m = 0; m <= l; ++m){
            poly.back().emplace_back(Poly(l, m));
            normCoef.back().emplace_back(sqrtq(((2*l+1)/(4*M_PI)) * fabsq(poly.back().back().getNeg())));
            fp >> aux;
            for(int i = 0; i < N; ++i){
                fp >> buf;
                poly.back().back().addVal(strtoflt128(buf.c_str(), NULL));
            }
        }
        fp.close();
    }*/
}

double SphereHarm::operator()(const int& l, const int& m, const double& tet, const double& phi) const
{
    if(std::abs(m) > l || l > lmax || l < 1) return 0;

    double aux = (normCoef[l][std::abs(m)]*poly[l][std::abs(m)](std::cos(tet)))*((m < 0) ? std::sin(-m*phi) : std::cos(m*phi));

    return aux;
}

double SphereHarm::gphi(const int& l, const int& m, const double& tet, const double& phi) const
{
    if(std::abs(m) > l || l > lmax || l < 1) return 0;

    return -m*normCoef[l][std::abs(m)]*poly[l][std::abs(m)](std::cos(tet))*((m < 0) ? std::cos(m*phi) : std::sin(m*phi));
}

double SphereHarm::gtet(const int& l, const int& m, const double& tet, const double& phi) const
{
    if(std::abs(m) > l || l > lmax || l < 1) return 0;
    else if(m == 0) return normCoef[l][0]*poly[l][1](std::cos(tet));
    else return normCoef[l][std::abs(m)] * (((std::abs(m)+1 <= l) ? poly[l][std::abs(m)+1](std::cos(tet)) : 0) + std::abs(m)*cotan(tet)*poly[l][std::abs(m)](std::cos(tet))) * ((m < 0) ? std::sin(-m*phi) : std::cos(m*phi));
}