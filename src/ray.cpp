#include "ray.hpp"

namespace ray {
    
    RayTraceOutput rayTrace(const Vector3D& x0, const Vector3D& v0, const RayTraceOpts& rtOpts)
    {
        Vector3D x = x0, v = v0;
        double tau = 0, r = rtOpts.dist(x), rmin;

        auto [n1, gradn] = rtOpts.gradN(x, r);
        rmin = r;

        do{
            // integrate position
            x += (rtOpts.dt*n1)*v;
            r = rtOpts.dist(x);
            if(r < rmin) rmin = r;

            // integrate velocity
            auto [n1, gradn] = rtOpts.gradN(x, r);
            v += (rtOpts.dt*sqr(n1))*(gradn - (v*gradn)*v);

            // integrate optical depth
            tau += (rtOpts.dt*n1)*rtOpts.optDepth(x, r);

            if(rtOpts.destroy(x, r)){
                tau = nb::inf;
                break;
            }

        }while(!rtOpts.exit(x, r));

        return {x, v, rmin, std::exp(-tau)};
    }

    RayTraceOutput rayTracePhase(const Vector3D& x0, const Vector3D& v0, const double& phi, const RayTraceOpts& rtOpts)
    {
        Vector3D x = x0, v = v0;
        double tau = 0, r = rtOpts.dist(x), rmin;
        int i = 0;

        auto [n1, gradn] = rtOpts.gradN(x, r);
        rmin = r;

        do{
            ++i;
            // integrate position
            x += (rtOpts.dt*n1)*v;
            r = rtOpts.dist(x);
            if(r < rmin) rmin = r;

            // integrate velocity
            auto [n1, gradn] = rtOpts.gradN(x, r);
            v += (rtOpts.dt*sqr(n1))*(gradn - (v*gradn)*v);

            // integrate optical depth
            tau += (rtOpts.dt*n1)*rtOpts.optDepth(x, r);

            if(rtOpts.destroy(x, r)){
                tau = nb::inf;
                break;
            }

        }while(!rtOpts.exit(x, r));

        return {x, v, rmin, std::exp(-tau), phi + i*rtOpts.dp};
    }
    
    
    
    
    
    
    
    Vector2D rayZplane(const Vector3D& x0, const Vector3D& v0, const double& z)
    {
        
        // check if ray is parallel to the plane
        if(v0[Z] == 0){
            #ifndef NO_CHECKS
                std::cout << "Warning: Ray is normal to the z-axis, no plane intersection is possible\n";
            #endif
            return (Vector2D){NAN,NAN};
        }

        // get time at which ray intersects with the plane
        double t = (z - x0[Z])/v0[Z];
        Vector2D res;

        #ifndef NO_CHECKS
        if(t < 0){
            std::cout << "Warning: Plane intersection lies behind the ray\n";
        }
        #endif

        // project the ray to the XY plane Z = z
        res[X] = x0[X] + v0[X]*t;
        res[Y] = x0[Y] + v0[Y]*t;

        return res;
    }

    Vector3D ray2sphere(const Vector3D& x0, const Vector3D& v0, const Vector3D& Rbb, const Vector3D& cc, const bool& in)
    {
        double c1, c2, c3, det;
        Vector3D dx = x0 - cc;

        // calculate coefficients of the quadratic polynomial
        c1 = sqr(v0[X]) + sqr(Rbb[Y]*v0[Y]) + sqr(Rbb[Z]*v0[Z]);
        c2 = v0[X]*dx[X] + sqr(Rbb[Y])*v0[Y]*dx[Y] + sqr(Rbb[Z])*v0[Z]*dx[Z];
        c3 = sqr(dx[X]) + sqr(Rbb[Y]*dx[Y]) + sqr(Rbb[Z]*dx[Z]) - sqr(Rbb[X]);

        // check if the solutions are real
        if((det = sqr(c2) - c1*c3) < 0){
            #ifndef NO_CHECK
                std::cout << (Vector3D){c1, c2, c3} << " " << det << std::endl;
                std::cout << "Warning: Ray doesn't intersect with the sphere\n";
            #endif
            return (Vector3D){NAN,NAN,NAN};
        }

        // return ray projection
        return x0 + v0*((-c2 + ((in) ? -std::sqrt(det) : std::sqrt(det)))/c1);
    }

    Vector3D ray2spherePhase(const Vector3D& x0, const Vector3D& v0, double& phi, const double& l, const Vector3D& Rbb, const Vector3D& cc, const bool& in)
    {
        double c1, c2, c3, det, t;
        Vector3D dx = x0 - cc;

        // calculate coefficients of the quadratic polynomial
        c1 = sqr(v0[X]) + sqr(Rbb[Y]*v0[Y]) + sqr(Rbb[Z]*v0[Z]);
        c2 = v0[X]*dx[X] + sqr(Rbb[Y])*v0[Y]*dx[Y] + sqr(Rbb[Z])*v0[Z]*dx[Z];
        c3 = sqr(dx[X]) + sqr(Rbb[Y]*dx[Y]) + sqr(Rbb[Z]*dx[Z]) - sqr(Rbb[X]);

        // check if the solutions are real
        if((det = sqr(c2) - c1*c3) < 0){
            #ifndef NO_CHECK
                std::cout << (Vector3D){c1, c2, c3} << " " << det << std::endl;
                std::cout << "Warning: Ray doesn't intersect with the sphere\n";
            #endif
            return (Vector3D){NAN,NAN,NAN};
        }

        t = (-c2 + ((in) ? -std::sqrt(det) : std::sqrt(det)))/c1;
        phi += 2*nb::pi/l;

        return x0 + v0*t;
    }
};