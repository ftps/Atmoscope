#ifndef PLANET_HPP
#define PLANET_HPP

#include <memory>

#include "utils.hpp"
#include "turbmap.hpp"

namespace planet {
    class Planet {
        public:
            Planet(const std::string& filename);
            
            void printParameters() const;

            // index of refraction functions
            refracOutput getRefrac(const Vector3D& pos, const double& r) const;
            refracOutput getRefracTurb(const Vector3D& pos, const double& r) const;

            // transmittance function
            double getOpticalDepth(const double& r) const;

            // extra functions for the ray tracer
            double barDistance(const Vector3D& pos) const;
            bool exitAtmosphere(const double& r) const;
            bool destroySurface(const double& r) const;


            std::string filename;
            double R, H, R_atm, n0;
            Vector3D beta, betaSqr;
            Matrix3D rtt, rti, dPhi;
            double latP, lonP;
        
        protected:
            std::unique_ptr<turb::TurbMap> turbMap;
            bool turbAvailable;
    };
};

#endif