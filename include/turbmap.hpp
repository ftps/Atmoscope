#ifndef TURBMAP_HPP
#define TURBMAP_HPP

#include "utils.hpp"

namespace turb {
    class TurbMap {
        public:
            TurbMap(const std::string& filename);
            Vector3D operator()(const double& tet, const double& phi) const;
            void readTurbFile(const std::string& filename);
            std::string getHeader();
        private:
            double maxTet, dt;
            uint imax, jmax;
            Vmap<double> map, gradTet, gradPhi;
            std::string headerFile;
    };
};


#endif