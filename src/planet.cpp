#include "planet.hpp"

namespace planet {
    Planet::Planet(const std::string& filename)
    {
        std::ifstream fp(filename);
        double ay, az;

        fp >> R;
        fp >> ay;
        fp >> az;

        beta = {1.0, 1.0/(1.0+ay), 1.0/(1.0+az)};
        betaSqr = beta^beta;

        fp >> H;
        fp >> n0;
        fp >> turbAvailable;

        if(turbAvailable){
            std::string turbfile;
            fp >> turbfile;
            turbMap = std::make_unique<turb::TurbMap>(turb::TurbMap((turbfile)));
        }

        R_atm = R + 20*H;

        this->filename = filename;
    }

    void Planet::printParameters() const
    {
        std::cout << filename << std::endl;
        std::cout << "\nR: " << R << " km\n";
        std::cout << "H: " << H << " km\n";
        std::cout << "n0: " << n0 << '\n';
        std::cout << "b: " << beta << "\n\n";
    }


    refracOutput Planet::getRefrac(const Vector3D& pos, const double& r) const
    {
        double aux = n0*std::exp((R-r)/H);
        return {1.0/(1.0 + aux), (betaSqr^pos)*(-aux/(H*r))};
    }

    /*Vector3D Planet::getRefracGrad(const Vector3D& pos, const double& r) const
    {
        double aux = -(n0/(H*r))*std::exp(-(r-R)/H);
        return (betaSqr^pos)*aux;
    }*/

    refracOutput Planet::getRefracTurb(const Vector3D& pos, const double& r) const
    {
        Vector3D newPos, turbGrad;
        Matrix3D pullBack;
        double tet, phi, rShort, n1 = n0*exp((R-r)/H);

        newPos = {pos[X]*std::cos(latP) - pos[Z]*beta[Z]*std::sin(latP), beta[Y]*pos[Y], pos[Z]*beta[Z]*std::cos(latP) + pos[X]*std::sin(latP)};
        rShort = std::sqrt(sqr(newPos[X]) + sqr(newPos[Y]));
        tet = std::atan2(newPos[Z], rShort);
        phi = std::atan2(newPos[Y], newPos[X]);
        if(phi < 0) phi += 2*nb::pi;

        turbGrad = turbMap->operator()(tet, phi)*n1;
        n1 = 1/(1+turbGrad[X]);
        turbGrad[X] = -turbGrad[X]/H;

        pullBack[X] = newPos/r;
        pullBack[Y] = {newPos[X]*newPos[Z]/(sqr(r)*rShort), newPos[X]*newPos[Y]/(sqr(r)*rShort), -rShort*beta[Z]/sqr(r)};
        pullBack[Z] = {-newPos[Y]/sqr(rShort), newPos[X]/sqr(rShort), 0};

        turbGrad = dPhi*(turbGrad*pullBack);

        return {n1, turbGrad};
    }






    double Planet::barDistance(const Vector3D& pos) const
    {
        return normP2(beta^pos);
    }

    bool Planet::exitAtmosphere(const double& r) const
    {
        return r > R_atm;
    }

    bool Planet::destroySurface(const double& r) const
    {
        return r < R;
    }
};