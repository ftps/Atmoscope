#ifndef RAY_HPP
#define RAY_HPP

#include "utils.hpp"

namespace ray {

    /*
        Structs for ray tracing
    */

    struct RayTraceOpts {
        double dt = 10;
        double l = 0;
        double dp = 0;
        double phi = 0;
        funcScaVec dist = [](const Vector3D& pos){ return std::sqrt(pos*pos); };
        funcScaVecSca optDepth = [](const Vector3D& pos, const double& r){ return 0; };
        funcRefrac gradN;
        funcExit exit, destroy;
    };

    struct RayTraceOutput {
        Vector3D x;
        Vector3D v;
        double rmin;
        double T = 1;
        double phi = 0;
    };



    /*
        Ray tracing function
    */

    RayTraceOutput rayTrace(const Vector3D& x0, const Vector3D& v0, const RayTraceOpts& rtOpts);

    /*
        Ray casting functions
    */

    // function that casts a ray to a plane normal to the z-axis
    Vector2D rayZplane(const Vector3D& x0, const Vector3D& v0, const double& z);

    // function that casts a ray to a spheroid, either the entering or exiting point
    // Rbb contains the x-radius at Rbb[X] and the beta values at Rbb[Y] and Rbb[Z]
    Vector3D ray2sphere(const Vector3D& x0, const Vector3D& v0, const Vector3D& Rbb, const Vector3D& cc = {0,0,0}, const bool& in = true);

    Vector3D ray2spherePH(const Vector3D& x0, const Vector3D& v0, double& phi, const Vector3D& Rbb, const Vector3D& cc, const double& l, const bool& in = true);
};


#endif