#ifndef FLASHMAP_HPP
#define FLASHMAP_HPP

#include "gnuplot-iostream.hpp"
#include "ray.hpp"
#include "planet.hpp"

#include <thread>
#include <chrono>

#define NTHREAD std::thread::hardware_concurrency()
#define R_std 6357.0

namespace flash {
    /*
        Struct for flash map options
    */

    struct FlashOpts {
        uint N = 250;
        double S = 50;
        double L = 1e8;
        double lat = 0;
        double lon = 0;
        bool turb = false;
        double dr_n = 1;
        double dphi_n = 1;
        double dt_n = 1;
        int multi = 0;
        ray::RayTraceOpts rtOpts;
    };

    struct FlashOutput {
        int i, j;
        double T;
        double phi;
    };

    /*
        Flash map classes
    */

    template <typename T>
    class fmap : public FlashOpts, public Vmap<T> {
        public:
            fmap();
            fmap(const FlashOpts& fOpts);
            void operator+=(const fmap& other);
            void multiAdd(const fmap& other);
            void normalize(double norm);
            double maximumVal() const;
            double maximumDist() const;

            double A;
            uint n;
            uint ntot;
    };

    class FlashMap : public planet::Planet, public fmap<double> {
        public:
            FlashMap(const std::string& filename, const FlashOpts& fOpts);
            ~FlashMap();

            void printFlashpars() const;
            void finalData(const std::chrono::duration<double>& interval) const;
            void plotMap(const std::string& filename = "", const bool& logScale = false) const;

            void setPosition(const double& lon, const double& lat);
            void setTurbulence(const bool& turb);

            void operator()();

            void singleMap();
            void multiMap();

            FlashOutput propagateRay(const Vector3D& x0s, const Vector3D& v0s) const;
            void map2detector(const Vector2D& xy, int& i, int& j) const;

            Vector2D getStartingRadius();
            inline int squareIntersect(const int& i1, const int& j1, const int& i2, const int& j2) const;
            inline void singleLoop(fmap& m, const double& cp, const double& sp, const double& rs, const bool& down) const;
            void mapThread(fmap& m, const double& phi_s, const double& phi_e, const Vector2D& rs, const int& iThread) const;

        private:
            double dphi, dr;
            const Vector3D v0 = Vector3D({0,0,1});
    };
};

#endif