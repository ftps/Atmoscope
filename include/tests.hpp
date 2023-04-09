#ifndef TESTS_HPP
#define TESTS_HPP

#include "gnuplot-iostream.hpp"
#include "flashmap.hpp"


double centralDifference(const funcScaSca& f, const double& x, const double& h);
double newtonRaphson(const funcScaSca& f, const double& x0, const double& a = 1, const double& e_rel = 1e-6, const double& e_abs = 1e-6, const uint& i_max = 500, const double& h = 1e-3);

namespace tests {
    void testUtils();
    void testCasts();
    void testRayTrace();
    void testPlanet();
    void testFlash();
    void testOblate();

    void runAll();


    void testFlashMap();
};

namespace geoOptics {
    void sphericalImages();
    void oblateImages(const double& lat = 0);
    void detectorSize();
    void detectorDistance(const bool& sph = false);
    void sphericalLaw();
    void oblateLaw();

    void runAll();
}

namespace turbAtmos {
    void Cn2(const double& n0 = 0.000278, const double& M = 0.0209391, const double& l0 = 100);
    void bendingAngle();
    void sphericalCase();
    void oblateCase();
    void prettyTitan();

    void runAll();
};

namespace waveOptics {
    void airyRing();
    void sphericalCase();
    void oblateCase();
    void prettyTitan();
    void horFluc();

    void runAll();
};

#endif
