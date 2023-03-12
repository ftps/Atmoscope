#ifndef TURB_HPP
#define TURB_HPP

#include "sphere.hpp"
#include <random>
#include <thread>

#define G53 0.9027452929509336112969
#define G16 5.56631600178023520425

#ifndef nthread
    #define nthread std::thread::hardware_concurrency()
#endif

struct turbOpts {
    //double n0 = 0.000277;
    double M = 0.0309391487;
    double R = 6357;
    double l0 = 250;
};

class Turbulence : public SphereHarm {
    public:
        Turbulence(const turbOpts& topts = turbOpts());
        void genCoef(const int& lplus = 50, const int& N = 100);
        void genMap(const int& Nphi, const double& tet_max = M_PI/8.0);
        void printMap(const std::string& filename);
        void plotMap(const std::string& filename = "", const int& div = 10);
        void genThread(const int& is, const int& ie, const int& Ntet, const double& dp);

    //private:
        bool isGen;
        double Co2, max_tet;
        int lmax_e;
        std::vector<std::vector<std::array<double,2>>> Clm;
        std::vector<std::vector<double>> map, g_tet, g_phi;
        VpairXY E, Er;
};

#endif