#include "test.hpp"

int main(int argc, char* argv[])
{
    turbOpts tEarth, tTitan;

    tTitan.M = 10/std::sqrt(1.4*8.3145*186/0.280135);;
    tTitan.R = 2516.375;
    tTitan.l0 = 600;

    std::vector<turbOpts> tOpts = {tEarth, tTitan};
    std::vector<std::string> planets = {"Earth", "Titan"};

    for(uint j = 0; j < tOpts.size(); ++j){
        Turbulence t(tOpts[j]);
        for(uint i = 1; i <= 30; ++i){
            std::cout << "Running atmosphere " << i << " for planet " << planets[j] << std::endl;
            t.genCoef(200, 200);
            t.genMap(500, M_PI/5.0);
            t.printMap("Data/" + planets[j] + "/atmos_" + std::to_string(i) + ".dat");
            t.plotMap("Data/" + planets[j] + "/atmos_" + std::to_string(i) + ".png", 12);
        }
    }

    return 0;
}