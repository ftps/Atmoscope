#include "turbmap.hpp"

namespace turb {
    TurbMap::TurbMap(const std::string& filename)
    {
        readTurbFile(filename);
    }

    void TurbMap::readTurbFile(const std::string& filename)
    {
        std::ifstream fp(filename);
        double aux;

        if(fp.bad()){
            exit(-2);
        }

        fp >> maxTet;
        fp >> imax;
        fp >> jmax;

        dt = 2*M_PI/(imax-1);

        map.clear();
        gradTet.clear();
        gradPhi.clear();

        for(uint i = 0; i < imax; ++i){
            map.push_back(std::vector<double>(0));
            gradTet.push_back(std::vector<double>(0));
            gradPhi.push_back(std::vector<double>(0));
            for(uint j = 0; j < jmax; ++j){
                fp >> aux;
                map.back().push_back(aux);
                fp >> aux;
                gradTet.back().push_back(aux);
                fp >> aux;
                gradPhi.back().push_back(aux);
            }
        }

        headerFile = filename.substr(0, filename.rfind('_'));
    }

    Vector3D TurbMap::operator()(const double& tet, const double& phi) const
    {
        Vector3D res;
        uint i, j;
        double tetN, phiN;

        if(std::abs(tet) > maxTet) return {1,0,0};
    
        i = std::floor((imax-1)*phi/(2*nb::pi));
        j = std::floor((tet/maxTet + 1.0)*0.5*(jmax-1));

        if(i == imax-1) --i;
        //if(j == jmax-1) --j;

        tetN = (tet + maxTet - j*dt)/dt;
        phiN = (phi - i*dt)/dt;

        if(i >= imax || j >= jmax || i < 0 || j < 0){
            std::cout << tet << ", " << phi << ", " << maxTet << std::endl;
            std::cout << i << ", " << j << std::endl;
            std::cout << imax << ", " << jmax << std::endl;
            std::cout << map.size() << ",  " << map[0].size() << std::endl;
        }

        res[X] = std::lerp(std::lerp(map[i][j], map[i+1][j], phiN),
                    std::lerp(map[i][j+1], map[i+1][j+1], phiN), tetN);
        res[Y] = std::lerp(std::lerp(gradTet[i][j], gradTet[i+1][j], phiN),
                    std::lerp(gradTet[i][j+1], gradTet[i+1][j+1], phiN), tetN);
        res[Z] = std::lerp(std::lerp(gradPhi[i][j], gradPhi[i+1][j], phiN),
                    std::lerp(gradPhi[i][j+1], gradPhi[i+1][j+1], phiN), tetN);

        return res;   
    }

    std::string TurbMap::getHeader()
    {
        return headerFile;
    }
};