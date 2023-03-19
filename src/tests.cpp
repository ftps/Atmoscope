#include "tests.hpp"

double centralDifference(const funcScaSca& f, const double& x, const double& h)
{
    return (f(x + h) - f(x - h))/(2*h);
}

double newtonRaphson(const funcScaSca& f, const double& x0, const double& a, const double& e_rel, const double& e_abs, const uint& i_max, const double& h)
{
    double x = x0, x_old = x0, aux_err;
    
    for(uint i = 0; i < i_max; ++i){
        x = x - a*(f(x)/centralDifference(f, x, h));
        aux_err = std::abs(x - x_old);
        if(std::abs(aux_err/x_old) < e_rel || aux_err < e_abs){
            return x;
        }
        x_old = x;
    }

    std::cout << "WARNING: Max iterations reached, last iteration will be returned.\n";
    return x; 
}


namespace tests
{
    void empty()
    {
        std::cout << "\n\n\n\n";
    }

    void testUtils()
    {
        Vector3D a = {2, 3, 4}, b = {0, 1, 1}, c = {1, 2, 3};
        Matrix3D M;
        std::vector<double> v = {1, 1.1, 1.2, 0.7, 0.9, 0.86, 1.12, 1.01, 0.88, 0.99};

        M[X] = a;
        M[Y] = b;
        M[Z] = c;

        empty();
        std::cout << "Testing Utils" << std::endl;
        std::cout << a << " " << b << " " << c << std::endl;
        std::cout << a + b << " " << c - b << std::endl;
        std::cout << 3.5 * a << " " << c / 2.3 << std::endl;
        std::cout << cross({1, 0, 0}, {0, 1, 0}) << " " << cross({0, 0, 1}, {0, 1, 0}) << std::endl << std::endl;

        std::cout << M << std::endl;
        std::cout << M * M << std::endl;
        std::cout << M * a << " " << b * M << std::endl << std::endl;

        std::cout << "1.141 " << sqr(1.141) << std::endl;
        std::cout << avg(v) << " " << var(v) << std::endl;
    }

    void testCasts()
    {
        Vector3D x0 = {2.3,1.4,-5}, v0 = {1.2, 0.3, -0.4};
        double z = 15;
        Vector3D Rbb = {15, 1, 1}, cc = {-10, 0, 0};

        empty();
        std::cout << "Testing Ray Casts" << std::endl;
        std::cout << ray::rayZplane(x0, {0,0,1}, z) << " " << ray::rayZplane(x0, v0, z) << std::endl;
        std::cout << ray::ray2sphere({0,0,0}, {1,0,0}, Rbb, cc) << " " << ray::ray2sphere({0,0,0}, {1,0,0}, Rbb, cc, false) << " " << ray::ray2sphere({0,0,20}, {1,0,0}, Rbb, cc) << std::endl;
    }

    void testRayTrace()
    {
        ray::RayTraceOpts rtOpts;
        double H = 8, R = 6357, n0 = 0.000277;
        Vector3D v0 = {0,0,1};
        funcScaSca bend = [&](const double& h){
            return n0*std::sqrt(2*nb::pi*R/H)*std::exp(-h/H);
        };
        VpairXY xy_calc, xy_numeric;
        Gnuplot gp;
        auto plot = gp.plotGroup();

        rtOpts.gradN = [&](const Vector3D& pos, const double& r){
            double aux = n0*std::exp(-(r-R)/H);
            return refracOutput({1.0/(1.0 + aux), pos*(-aux/(H*r))});
        };

        rtOpts.exit = [&](const Vector3D& pos, const double& r){
            return r > (R + 15*H);
        };

        rtOpts.destroy = [&](const Vector3D& pos, const double& r){
            return r < R;
        };

        rtOpts.dt = 1;

        empty();
        std::cout << "Testing ray trace" << std::endl;
        for(double i = R + 14*H; i > R; i -= 1){
            Vector3D x0 = ray::ray2sphere({i, 0, -5*R}, v0, {R+15*H, 1, 1});

            auto [x, v, rmin, T, phi] = ray::rayTrace(x0, v0, rtOpts);

            if(T == 0){
                std::cout << "Contact" << std::endl;
                break;
            }

            xy_calc.emplace_back(rmin - R, bend(rmin-R)*rad2deg);
            xy_numeric.emplace_back(rmin - R, std::abs(atan2(v[X], v[Z]))*rad2deg);
        }

        plot.add_plot1d(xy_calc, "with linespoints title 'Approximation' lw 3 ps 2 pt 7");
        plot.add_plot1d(xy_numeric, "with linespoints title 'Numeric' lw 3 ps 2 pt 7");

        gp << "set terminal png size 1600,1200 enhanced font 'Verdana,30' enhanced\nset encoding utf8\n";
        gp << "set termoption enhanced\n";
        gp << "set output 'Results/testRayTrace.png\n";
        gp << "set logscale y\n";
        gp << "set xlabel 'h_{min} - km'\n";
        gp << "set ylabel 'Δθ - º' enhanced\n";
        gp << plot;

    }

    void testPlanet()
    {
        ray::RayTraceOpts rtOpts;
        planet::Planet p("Planets/EarthSpherical.dat");
        Vector3D v0 = {0,0,1};
        funcScaSca bend = [&](const double& h){
            return p.n0*std::sqrt(2*nb::pi*p.R/p.H)*std::exp(-h/p.H);
        };
        VpairXY xy_calc, xy_numeric;
        Gnuplot gp;
        auto plot = gp.plotGroup();

        rtOpts.gradN = [&](const Vector3D& pos, const double& r){
            return p.getRefrac(pos, r);
        };

        rtOpts.exit = [&](const Vector3D& pos, const double& r){
            return p.exitAtmosphere(r);
        };

        rtOpts.destroy = [&](const Vector3D& pos, const double& r){
            return p.destroySurface(r);
        };

        rtOpts.dist = [&](const Vector3D& pos){
            return p.barDistance(pos);
        };

        rtOpts.dt = 1;

        empty();
        std::cout << "Testing planet" << std::endl;
        for(double i = p.R + 14*p.H; i > p.R; i -= 1){
            Vector3D x0 = ray::ray2sphere({i, 0, -5*p.R}, v0, {p.R+15*p.H, 1, 1});

            auto [x, v, rmin, T, phi] = ray::rayTrace(x0, v0, rtOpts);

            if(T == 0){
                std::cout << "Contact" << std::endl;
                break;
            }

            xy_calc.emplace_back(rmin - p.R, bend(rmin-p.R)*rad2deg);
            xy_numeric.emplace_back(rmin - p.R, std::abs(atan2(v[X], v[Z]))*rad2deg);
        }

        plot.add_plot1d(xy_calc, "with linespoints title 'Approximation' lw 3 ps 2 pt 7");
        plot.add_plot1d(xy_numeric, "with linespoints title 'Numeric' lw 3 ps 2 pt 7");

        gp << "set terminal png size 1600,1200 enhanced font 'Verdana,30' enhanced\nset encoding utf8\n";
        gp << "set output 'Results/testPlanet.png\n";
        gp << "set termoption enhanced\n";
        gp << "set logscale y\n";
        gp << "set xlabel 'h_{min} - km' enhanced\n";
        gp << "set ylabel 'Δθ - º' enhanced\n";
        gp << plot;
    }

    void testFlash()
    {
        flash::FlashOpts fopts;
        flash::FlashMap f("Planets/EarthSpherical.dat", fopts);

        empty();
        std::cout << "Testing Flash Map\n" << std::endl;

        f();

        std::cout << "\nNumerical: " << f.maximumVal() << ", Theoretical: " << 7.0509987*f.H*f.N/f.S << std::endl;

        f.plotMap("Results/testFlash.png");
    }

    void testOblate()
    {
        flash::FlashOpts fopts;
        fopts.S = 100;
        flash::FlashMap f("Planets/Earth.dat", fopts);

        empty();
        std::cout << "Testing Oblate Case\n" << std::endl;

        f();

        std::cout << "\nTheoretical Spherical: " << 7.0509987*f.H*f.N/f.S << ", Oblate Case: " << f.maximumVal() << std::endl;

        f.plotMap("Results/testOblate1.png");

        f.setPosition(0, 60);
        f();

        std::cout << "\nTheoretical Spherical: " << 7.0509987*f.H*f.N/f.S << ", Oblate Case: " << f.maximumVal() << std::endl;

        f.plotMap("Results/testOblate2.png");
    }




    void runAll()
    {
        testUtils();
        testCasts();
        testRayTrace();
        testPlanet();
        testFlash();
        testOblate();
    }
};


namespace geoOptics {
    void sphericalImages()
    {
        flash::FlashOpts fOpts;
        std::vector<std::string> planets = {"Earth", "Titan"};
        std::ofstream fp("Results/DataFiles/geoOptics_Amp.dat");
        
        fOpts.S = 1.51;
        fOpts.N = 151;
        //fOpts.lat = 90;

        for(const std::string& p : planets){
            flash::FlashMap f("Planets/" + p + ".dat", fOpts);

            f();

            f.plotMap("Results/geoOptics_" + p + "_spherical.png");
            fp << "Planet: " << p << ", Lat: " << fOpts.lat << "º" << std::endl;
            fp << "Amplification: " << f.maximumVal() << ", Detector Size: " << f.S/f.N << " km, Diamond size: " << 2*f.maximumDist() << " km" << std::endl << std::endl;
        }
    }

    void oblateImages(const double& lat)
    {
        flash::FlashOpts fOpts;
        std::vector<std::string> planets = {"Earth", "Titan"};
        std::ofstream fp("Results/DataFiles/geoOptics_Amp.dat", std::ios_base::app);
        
        fOpts.lat = lat;
        fOpts.S = 20;
        fOpts.N = 750;

        for(const std::string& p : planets){
            flash::FlashMap f("Planets/" + p + ".dat", fOpts);

            f();

            f.plotMap("Results/geoOptics_" + p + "_obl_" + std::to_string(lat) + ".png");
            fp << "Planet: " << p << ", Lat: " << fOpts.lat << "º" << std::endl;
            fp << "Amplification: " << f.maximumVal() << ", Detector Size: " << f.S/f.N << " km, Diamond size: " << 2*f.maximumDist() << " km" << std::endl << std::endl;
        }
    }

    void detectorSize()
    {
        Gnuplot gp;
        auto plot = gp.plotGroup();
        VpairXY xy;
        flash::FlashOpts fOpts;
        std::vector<std::string> planets = {"Earth", "Titan", "Mars", "Uranus", "Saturn"};
        //std::vector<double> k0Vals;
        //std::ofstream fp("Results/DataFiles/sphConstant.dat");
        double dS = std::pow(10, -1.0/5.0), S0 = 15.1;//, k0ci;

        //fOpts.S = 1.51;
        fOpts.N = 151;
        fOpts.lat = 90;

        for(const std::string& p : planets){
            flash::FlashMap f("Planets/" + p + ".dat", fOpts);
            xy.clear();
            
            for(double S = S0; S > 0.1; S *= dS){
                f.S = S;
                f();
                xy.emplace_back(1000.0*f.S/f.N, f.maximumVal());

                //k0Vals.push_back(f.maximumVal()*f.beta[Y]*f.S/(f.N*f.H));

                //fp << "H: " << f.H << " km, d0: " << f.S/f.N << " m" << std::endl;
                //fp << "k0: " << k0Vals.back() << std::endl;
            }

            plot.add_plot1d(xy, "with linespoints title '" + p + "' lw 3 ps 2 pt 7");
        }

        gp << "set terminal png size 1600,1200 enhanced font 'Verdana,30' enhanced\n";
        gp << "set termoption enhanced\n";
        gp << "set output 'Results/geoOptics_detSize.png\n";
        gp << "set logscale y\nset logscale x\n";
        gp << "set xlabel 'd_{0} - m'\n";
        gp << "set ylabel 'Amplification'\n";
        gp << plot;

        //k0ci = stdDev(k0Vals)*tStudentCritical95(k0Vals.size())/std::sqrt(k0Vals.size());

        //fp << "k0 average: " << avg(k0Vals) << ", k1 std: " << stdDev(k0Vals) << ", c.i.: " << k0ci << std::endl;
    }

    void detectorDistance(const bool& sph)
    {
        Gnuplot gp;
        auto plot = gp.plotGroup();
        VpairXY xy;
        flash::FlashOpts fOpts;
        std::vector<std::string> planets = {"Earth", "Titan", "Mars", "Uranus", "Saturn"};
        double dL = std::pow(10, 1.0/5.0), L0 = 1e7;

        fOpts.S = 0.151;
        fOpts.N = 151;
        if(!sph) fOpts.lat = 90;

        for(const std::string& p : planets){
            flash::FlashMap f("Planets/" + p + ((sph) ? "Spherical.dat" : ".dat"), fOpts);
            xy.clear();
            
            for(double L = L0; L < 1e11; L *= dL){
                f.L = L;
                f();
                xy.emplace_back(f.L, f.maximumVal());
            }

            plot.add_plot1d(xy, "with linespoints title '" + p + "' lw 3 ps 2 pt 7");
        }

        gp << "set terminal png size 1600,1200 enhanced font 'Verdana,30' enhanced\n";
        gp << "set termoption enhanced\n";
        gp << "set output 'Results/geoOptics_detDist" << ((sph) ? "Spherical.png'\n" : ".png'\n");
        gp << "set logscale y\nset logscale x\n";
        gp << "set xlabel 'L - km'\n";
        gp << "set ylabel 'Amplification'\n";
        gp << plot;
    }

    void sphericalLaw()
    {
        flash::FlashOpts fOpts;
        std::vector<std::string> planets = {"Earth", "Titan", "Mars", "Uranus", "Saturn"};
        std::vector<double> k0Vals;
        std::ofstream fp("Results/DataFiles/sphConstant.dat");
        double dS = std::pow(10, -1.0/5.0), S0 = 15.1, k0ci;

        //fOpts.S = 1.51;
        fOpts.N = 151;
        fOpts.lat = 90;

        for(const std::string& p : planets){
            flash::FlashMap f("Planets/" + p + "Spherical.dat", fOpts);
            
            for(double S = S0; S > 0.1; S *= dS){
                f.S = S;
                f();
                //xy.emplace_back(1000.0*f.S/f.N, f.maximumVal());

                k0Vals.push_back(f.maximumVal()*f.beta[Y]*f.S/(f.N*f.H));

                fp << "H: " << f.H << " km, d0: " << f.S/f.N << " m" << std::endl;
                fp << "k0: " << k0Vals.back() << std::endl;
            }
        }

        k0ci = stdDev(k0Vals)*tStudentCritical95(k0Vals.size())/std::sqrt(k0Vals.size());

        fp << "k0 average: " << avg(k0Vals) << ", k1 std: " << stdDev(k0Vals) << ", c.i.: " << k0ci << std::endl;
    }

    void oblateLaw()
    {
        Gnuplot gp;
        auto plot = gp.plotGroup();
        VpairXY xy;
        std::vector<double> k1Vals, k2Vals;
        double dh = std::pow(10.0, 1.0/3.0), lat, k1ci, k2ci;
        flash::FlashOpts fOpts;
        std::ofstream fp("Results/DataFiles/oblConstant.dat");
        uint i = 0;

        fOpts.S = 3.001;
        fOpts.N = 151;

        flash::FlashMap f("Planets/Earth.dat", fOpts);
        double Rb = f.R/f.beta[Y];

        funcScaSca lt2t = [&](const double& theta){
            double aux = (theta == 90) ? 1/f.beta[Y] : (std::sqrt((1 + sqr(std::tan(theta*deg2rad)))/(1 + sqr(f.beta[Y]*std::tan(theta*deg2rad)))));

            return Rb*(1 - f.beta[Y]*aux);
        };

        for(double l = 0; l <= 1; l += 0.25){
            xy.clear();
            lat = (l == 0) ? 90 : newtonRaphson([&](const double& t){ return lt2t(t) - l; }, 80);
            f.setPosition(0, lat);
            for(double h = 1.0; h <= 10000.01; h *= dh){
                f.S = 4.5*lt2t(lat);
                if(f.S < 50*h/1000) f.S = 150*h/1000;
                f.N = std::floor(1000*f.S/h);
                f.S = ((f.N % 2 == 0) ? ++f.N : f.N)*h/1000.0;

                f();
                f.plotMap("Results/aux_" + std::to_string(i++) + ".png");

                xy.emplace_back(h, f.maximumVal());
                if(l != 0 && h < 200){
                    k1Vals.push_back(f.maximumVal()*std::pow(h/1000, 2.0/3.0)*l/f.H);
                    k2Vals.push_back(f.maximumVal()*std::pow(h/1000, 2.0/3.0)*std::pow(l, 1.0/3.0)/f.H);

                    fp << "H: " << f.H << " km, d0: " << h << " m, Δr: " << l << " km" << std::endl;
                    fp << "k1: " << k1Vals.back() << ", k2: " << k2Vals.back() << std::endl << std::endl;
                }
            }
            plot.add_plot1d(xy, "with linespoints title 'Δr = " + std::to_string(l) + " km' lw 3 ps 2 pt 7 enhanced");
        }

        gp << "set terminal png size 1600,1200 enhanced font 'Verdana,30' enhanced\nset encoding utf8\n";
        gp << "set termoption enhanced\n";
        gp << "set output 'Results/geoOptics_oblLaw.png'\n";
        gp << "set logscale y\nset logscale x\n";
        gp << "set xlabel 'd_0 - m'\n";
        gp << "set ylabel 'Amplification'\n";
        gp << plot;

        k1ci = stdDev(k1Vals)*tStudentCritical95(k1Vals.size())/std::sqrt(k1Vals.size());
        k2ci = stdDev(k2Vals)*tStudentCritical95(k2Vals.size())/std::sqrt(k2Vals.size());

        fp << "k1 average: " << avg(k1Vals) << ", k1 std: " << stdDev(k1Vals) << ", c.i.: " << k1ci << std::endl;
        fp << "k2 average: " << avg(k2Vals) << ", k2 std: " << stdDev(k2Vals) << ", c.i.: " << k2ci << std::endl;

    }






    void runAll()
    {
        std::vector<double> lats = {85, 75};

        sphericalImages();
        for(const double& lat : lats){
            oblateImages(lat);
        }

        detectorSize();
        detectorDistance();
        detectorDistance(true);
        sphericalLaw();
        oblateLaw();
    }
};




namespace turbAtmos {
    void Cn2(const double& l0)
    {
        Gnuplot gp;
        auto plot = gp.plotGroup();
        VpairXY xy1, xy2;
        double cn2 = 2.0*sqr(0.000278*sqr(0.0209391))*std::pow(2*nb::pi/l0, 2.0/3.0)/(5.0*0.90274529), H = 8.0;
        funcScaSca hTeo = [&](const double& h){
            return cn2*std::exp(-2*h/H);
        };
        funcScaSca hExp = [&](const double& h){
            if(h <= 2.13){
                return std::pow(10, -10.7025 - 4.3507*h + 0.8141*h*h);
            }
            else if(h <= 10.34){
                return std::pow(10, -16.2897 + 0.0335*h - 0.0134*h*h);
            }
            return std::pow(10, -17.0557 - 0.0449*h - 0.0005*h*h + 0.6181*std::exp(-0.5*sqr((h - 15.5617)/3.466)));
        };

        for(double h = 0; h < 100; h += 0.1){
            xy1.emplace_back(hTeo(h), h);
            xy2.emplace_back(hExp(h), h);
        }

        plot.add_plot1d(xy1, "with lines title 'Simple Model' lw 3");
        plot.add_plot1d(xy2, "with lines title 'Data Fit' lw 3");

        gp << "set terminal png size 1600,1000 enhanced font 'Verdana,30' enhanced\nset encoding utf8\n";
        gp << "set termoption enhanced\n";
        gp << "set output 'Results/turbAtmos_Cn2test.png'\n";
        gp << "set logscale x\n";
        gp << "set xlabel 'C_n^2 - m^{2/3}'\n";
        gp << "set ylabel 'h - km'\n";
        gp << plot;
    }

    void bendingAngle()
    {
        Gnuplot gp, gp2, gp3;
        auto plot = gp.plotGroup(), plot2 = gp2.plotGroup(), plot3 = gp3.plotGroup();
        VpairXY xy, xy_clear, xyx, xyy, xyx_c, xyy_c, xyrad, xyrad_c;
        planet::Planet p("Planets/EarthSpherical.dat");
        ray::RayTraceOpts rtOpts;
        Vector3D v0 = {0,0,1};
        double r = 6402.83, aux;

        p.rtt[X] = {1,0,0};
        p.rtt[Y] = {0,1,0};
        p.rtt[Z] = {0,0,1};
        p.rti = p.rtt;

        rtOpts.dist = [&](const Vector3D& pos){
            return p.barDistance(pos);
        };
        rtOpts.exit = [&](const Vector3D& pos, const double& r){
            return p.exitAtmosphere(r);
        };
        rtOpts.destroy = [&](const Vector3D& pos, const double& r) {
            return p.destroySurface(r);
        };

        for(double phi = 0; phi < 2*nb::pi+0.01; phi += 0.01){
            Vector3D x0 = {r*std::cos(phi), r*std::sin(phi), -p.R_atm};
            x0 = ray::ray2sphere(x0, v0, {p.R_atm, p.beta[Y], p.beta[Z]});
            
            rtOpts.gradN = [&](const Vector3D& pos, const double& r){
                return p.getRefracTurb(pos, r);
            };
            auto [x, v, auxr, auxT, auxp] = ray::rayTrace(x0, v0, rtOpts);
            xy.emplace_back(phi, 180*std::atan2(std::sqrt(sqr(v[X]) + sqr(v[Y])), v[Z])/nb::pi);

            Vector2D xy_plane = ray::rayZplane(x, v, 1e8);
            xyx.emplace_back(phi, xy_plane[X]);
            xyy.emplace_back(phi, xy_plane[Y]);
            aux = std::atan2(xy_plane[Y], xy_plane[X]);
            xyrad.emplace_back(phi, aux + ((aux < 0) ? 2*nb::pi : 0));

            rtOpts.gradN = [&](const Vector3D& pos, const double& r){
                return p.getRefrac(pos, r);
            };
            auto [x2, v2, auxr2, auxT2, auxp2] = ray::rayTrace(x0, v0, rtOpts);
            xy_clear.emplace_back(phi, 180*std::atan2(std::sqrt(sqr(v2[X]) + sqr(v2[Y])), v2[Z])/nb::pi);

            xy_plane = ray::rayZplane(x2, v2, 1e8);
            xyx_c.emplace_back(phi, xy_plane[X]);
            xyy_c.emplace_back(phi, xy_plane[Y]);
            aux = std::atan2(xy_plane[Y], xy_plane[X]);
            xyrad_c.emplace_back(phi, aux + ((aux < 0) ? 2*nb::pi : 0));
        }

        plot.add_plot1d(xy, "with lines title 'Turbulent' lw 3");
        plot.add_plot1d(xy_clear, "with lines title 'Ideal' lw 3");

        plot2.add_plot1d(xyx, "with lines title 'x - Turbulent' lw 3");
        plot2.add_plot1d(xyy, "with lines title 'y - Turbulent' lw 3");
        plot2.add_plot1d(xyx_c, "with lines title 'x - Ideal' lw 3");
        plot2.add_plot1d(xyy_c, "with lines title 'y - Ideal' lw 3");
        

        plot3.add_plot1d(xyrad, "with lines title 'Turbulent'");
        plot3.add_plot1d(xyrad_c, "with lines title 'Ideal'");

        gp << "set terminal png size 1600,1200 enhanced font 'Verdana,30' enhanced\nset encoding utf8\n";
        gp << "set termoption enhanced\n";
        gp << "set output 'Results/turbAtmos_bendAngle.png'\n";
        //gp << "set logscale y\n";
        gp << "set xtics (\"0\" 0, \"π/2\" pi/2, \"π\" pi, \"3π/2\" 3*pi/2, \"2π\" 2*pi)\n";
        gp << "set xlabel 'φ - rad'\n";
        gp << "set ylabel 'Δθ - º'\n";
        gp << "set grid x\n set grid y\n";
        gp << plot;

        gp2 << "set terminal png size 1600,1200 enhanced font 'Verdana,30' enhanced\nset encoding utf8\n";
        gp2 << "set termoption enhanced\n";
        gp2 << "set output 'Results/turbAtmos_rayHit.png'\n";
        //gp << "set logscale y\n";
        gp2 << "set xtics (\"0\" 0, \"π/2\" pi/2, \"π\" pi, \"3π/2\" 3*pi/2, \"2π\" 2*pi)\n";
        gp2 << "set xlabel 'φ - rad'\n";
        gp2 << "set ylabel 'x/y - km'\n";
        gp2 << "set grid x\n set grid y\n";
        gp2 << plot2;

        gp3 << "set terminal png size 1600,1200 enhanced font 'Verdana,30' enhanced\nset encoding utf8\n";
        gp3 << "set termoption enhanced\n";
        gp3 << "set output 'Results/turbAtmos_phiAngle.png'\n";
        //gp << "set logscale y\n";
        gp3 << "set xtics (\"0\" 0, \"π/2\" pi/2, \"π\" pi, \"3π/2\" 3*pi/2, \"2π\" 2*pi)\n";
        gp3 << "set ytics (\"0\" 0, \"π/2\" pi/2, \"π\" pi, \"3π/2\" 3*pi/2, \"2π\" 2*pi)\n";
        gp3 << "set xlabel 'φ - rad'\n";
        gp3 << "set ylabel \"φ' - rad\"\n";
        gp3 << "set grid x\n set grid y\n";
        gp3 << plot3;
    }

    void sphericalCase()
    {
        flash::FlashOpts fOptsE, fOptsT;
        fOptsE.S = 2.01;
        fOptsE.N = 201;
        fOptsE.turb = true;
        fOptsE.lat = 90;
        fOptsT.S = 40.2;
        fOptsT.N = 201;
        fOptsT.turb = true;
        fOptsT.lat = 90;
        std::vector<flash::FlashOpts> opts = {fOptsE, fOptsT};
        std::vector<std::string> pla = {"Earth", "Titan"};

        for(uint i = 0; i < pla.size(); ++i){
            flash::FlashMap f("Planets/" + pla[i] + ".dat", opts[i]);
            f();
            f.plotMap("Results/turbAtmos_" + pla[i] + "Single.png");
            f.multi = 30;
            f();
            f.plotMap("Results/turbAtmos_" + pla[i] + "Multi.png");

            if(i == 0){
                Gnuplot gp;
                auto plot = gp.plotGroup();
                VpairXY xy1, xy2;
                uint N = opts[i].N;
                double S = opts[i].S;
                double maxVal = f.maximumVal();

                for(uint k = 0; k < N; ++k){
                    xy1.emplace_back(k*S/(N-1) - S/2, f[(N-1)/2][k]);
                }

                f.turb = false;
                //f.N = 4*N;
                f();

                N = f.N;
                S = f.S;
                for(uint k = 0; k < N; ++k){
                    xy2.emplace_back((k+1)*S/N - S/2, f[(N-1)/2][k]);
                }

                plot.add_plot1d(xy1, "with lines title 'Turbulent' lw 3");
                plot.add_plot1d(xy2, "with lines title 'Ideal' lw 3");

                gp << "set terminal png size 1600,1000 enhanced font 'Verdana,30' enhanced\nset encoding utf8\n";
                gp << "set termoption enhanced\n";
                gp << "set output 'Results/turbAtmos_horBright.png'\n";
                gp << "set xlabel 'x - km'\n";
                gp << "set ylabel 'A'\n";
                gp << "set yrange [0:" << std::to_string(3*maxVal) << "]\n";
                gp << "set xrange [-1:1]\n";
                gp << plot;
            }
        }
    }

    void oblateCase()
    {
        flash::FlashOpts fOptsE, fOptsT;
        fOptsE.S = 16.08;
        fOptsE.N = 201;
        fOptsE.turb = true;
        fOptsE.lat = 75;
        fOptsT.S = 40.2;
        fOptsT.N = 201;
        fOptsT.turb = true;
        fOptsT.lat = 75;
        std::vector<flash::FlashOpts> opts = {fOptsE, fOptsT};
        std::vector<std::string> pla = {"Earth", "Titan"};

        for(uint i = 0; i < pla.size(); ++i){
            flash::FlashMap f("Planets/" + pla[i] + ".dat", opts[i]);
            f();
            f.plotMap("Results/turbAtmos_" + pla[i] + "OblSingle.png");
            f.multi = 30;
            f();
            f.plotMap("Results/turbAtmos_" + pla[i] + "OblMulti.png");
        }
    }

    void testFlashMap()
    {
        /*flash::FlashOpts fOpts;
        fOpts.S = 201;
        fOpts.N = 201;
        fOpts.turb = true;
        fOpts.multi = 0;
        fOpts.lat = 45;
        //fOpts.lat = 15;
        flash::FlashMap f("Planets/Titan.dat", fOpts);
        f();
        f.plotMap();*/

        flash::FlashOpts fOpts;
        fOpts.L = 1e8;
        fOpts.S = 50;
        fOpts.N = 501;
        fOpts.turb = false;
        fOpts.multi = 0;
        fOpts.lat = 90;
        fOpts.waveOptics = true;
        fOpts.rtOpts.l = 0.1;
        //fOpts.lat = 15;
        flash::FlashMap f("Planets/Earth.dat", fOpts);
        f();
        f.plotMap("", true);
    }
};