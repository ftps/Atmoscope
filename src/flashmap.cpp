#include "flashmap.hpp"

//#undef NTHREAD
//#define NTHREAD 1

namespace flash {

    template<typename T>
    fmap<T>::fmap() { }

    template<typename T>
    fmap<T>::fmap(const FlashOpts& fOpts) : FlashOpts(fOpts)
    {
        if(N % 2 == 0) ++N;
        this->resize(N, std::vector<T>(N, 0));
        
        A = 0;
        n = 0;
        ntot = 0;
    }

    template<typename T>
    void fmap<T>::operator+=(const fmap<T>& other)
    {
        for(luint i = 0; i < this->size(); ++i){
            for(luint j = 0; j < this->size(); ++j){
                this->at(i)[j] += other[i][j];
            }
        }

        A += other.A;
        n += other.n;
        ntot += other.ntot;
    }

    template<typename T>
    void fmap<T>::multiAdd(const fmap& other)
    {
        for(luint i = 0; i < this->size(); ++i){
            for(luint j = 0; j < this->size(); ++j){
                this->at(i)[j] += other[i][j] + other[j][this->size() - 1 - i] + other[this->size() - 1 - i][this->size() - 1 - j] + other[j][this->size() - 1 - i];
            }
        }

        A += other.A;
        n += other.n;
        ntot += other.ntot;
    }

    template<typename T>
    void fmap<T>::normalize(double norm)
    {
        for(luint i = 0; i < this->size(); ++i){
            for(luint j = 0; j < this->size(); ++j){
                this->at(i)[j] *= norm;
            }
        }
    }

    template<>
    double fmap<double>::maximumVal() const
    {
        double max = 0, aux;

        for(const std::vector<double>& v : *this){
            aux = *std::max_element(v.begin(), v.end());
            if(aux > max) max = aux;
        }

        return max;
    }

    template<>
    double fmap<double>::maximumDist() const
    {
        double max = 0;
        uint im = 0, jm = 0;

        for(uint i = 0; i < this->size(); ++i){
            for(uint j = 0; j < this->size(); ++j){
                if(this->at(i)[j] > max){
                    max = this->at(i)[j];
                    im = i;
                    jm = j;
                }
            }
        }

        return std::sqrt(sqr(im - (N-1)/2) + sqr(jm - (N-1)/2))*S/N;
    }


    FlashMap::FlashMap(const std::string& filename, const FlashOpts& fOpts) : planet::Planet(filename), fmap(fOpts)
    {
        // define function according to the planet's parameters
        rtOpts.dist = [&](const Vector3D& pos){
            return this->barDistance(pos);
        };
        rtOpts.exit = [&](const Vector3D& pos, const double& r){
            return this->exitAtmosphere(r);
        };
        rtOpts.destroy = [&](const Vector3D& pos, const double& r) {
            return this->destroySurface(r);
        };

        setTurbulence(turb);

        // wave optics sets
        rtOpts.dp = (rtOpts.l > 0) ? 2*nb::pi*std::fmod(rtOpts.dt/rtOpts.l, 1.0) : 0;

        // rotation matrices
        setPosition(lon, lat);
    }

    FlashMap::~FlashMap()
    {
        fs::path path{"data_dump.dat"};
        if(fs::exists(path)){
            fs::remove("data_dump.dat");
        }
    }

    void FlashMap::printFlashpars() const
    {
        std::cout << "Res: " << N << "x" << N << " pixels\n";
        std::cout << "Size: " << S << " km\n";
        std::cout << "Pixel size: " << 1000*S/N << " m\n";
        std::cout << "Distance: " << L << " km\n";
        std::cout << "Lat: " << lat << "º, Lon: " << lon << "º\n\n";
    }

    void FlashMap::finalData(const std::chrono::duration<double>& interval) const
    {
        std::cout << "Computational time: " << interval.count() << " s\n";
        std::cout << "Number of rays in the detector: " << n << std::endl;
        std::cout << "Number of total rays traced: " << ntot << std::endl;
        std::cout << "Hit percentage: " << (100.0*n)/ntot << "%\n";
        std::cout << "Maximum Amplification: " << maximumVal() << std::endl;
    }

    void FlashMap::plotMap(const std::string& filename, const bool& logScale) const
    {
        Gnuplot gp;
        std::ostringstream ss;
        std::ofstream fp("data_dump.dat");

        if(fp.bad()){
            std::cout << "Error: Cannot open auxiliar file for plotting the flash map\n";
            return;
        }

        for(const std::vector<double>& v : *this){
            for(const double& val : v){
                fp << val << " ";
            }
            fp << '\n';
        }
        fp.close();

        gp << "reset\nunset key\n";
        if(filename == ""){
            gp << "set terminal wxt size 800,600 enhanced font 'Verdana,15' persist\n";
        }
        else{
            gp << "set terminal png size 1600,1200 enhanced font 'Verdana,30'\n";
            gp << "set output '" << filename << "'\n";
        }

        gp << "set style line 11 lc rgb '#808080' lt 1\n";
        gp << "set border 3 front ls 11\n";
        gp << "set tics nomirror out scale 0.75\n";
        gp << "set palette defined (0 0 0 0, 1 1 0 0, 2 1 1 1)\n";
        gp << "set autoscale noextend\n";

        ss << -S/2 << ":" << S/2;
        gp << "set xrange [" + ss.str() + "]\nset xlabel 'x - km'\n";
        gp << "set yrange [" + ss.str() + "]\nset ylabel 'y - km'\n";

        ss.str("");
        ss.clear();
        ss << "+1)*" << S/N << "-" << S/2;
        gp << "set view map\nset size ratio -1\n";

        gp << "set pm3d interpolate 2,2\nsplot 'data_dump.dat' matrix using (($1" + ss.str() + "):(($2" + ss.str() + "):" + ((logScale) ? "(log10($3))" : "($3)") + " with pm3d\n";
    }


    void FlashMap::setPosition(const double& lon, const double& lat)
    {
        Matrix3D rpp;

        this->lon = lon;
        this->lat = lat;
        this->lonP = lon*deg2rad;
        this->latP = lat*deg2rad;

        // auxiliary matrix with longitudonal rotation
        rpp[X] = {1,0,0};
        rpp[Y] = {0,std::cos(lonP), -std::sin(lonP)};
        rpp[Z] = {0,std::sin(lonP), std::cos(lonP)};

        // rotation matrix with latitudonal rotation
        rtt[X] = {std::cos(latP),0,-std::sin(latP)};
        rtt[Y] = {0,1,0};
        rtt[Z] = {std::sin(latP),0,std::cos(latP)};

        // inverse rotation matrix with latitudonal rotation
        rti[X] = {std::cos(latP),0,std::sin(latP)};
        rti[Y] = {0,1,0};
        rti[Z] = {-std::sin(latP),0,std::cos(latP)};

        // full rotation of longitude + latitude
        rtt = rtt*rpp;
        // reversing longitudonal rotation
        rpp[Y][Z] *= -1;
        rpp[Z][Y] *= -1;
        // full inverse rotation of latitude + longitude
        rti = rpp*rti;

        // vector pullback
        dPhi[X] = {std::cos(latP), 0, std::sin(latP)};
        dPhi[Y] = {0, beta[Y], 0};
        dPhi[Z] = {-beta[Z]*std::sin(latP), 0, beta[Z]*std::cos(latP)};
    }

    void FlashMap::setTurbulence(const bool& turb)
    {
        this->turb = turb;
        if(turb && turbAvailable){
            rtOpts.gradN = [&](const Vector3D& pos, const double& r){
                return this->getRefracTurb(pos, r);
            };
        }
        else{
            if(turb){
                std::cout << "\nWARNING: No turbulent atmosphere available for the current planet, an ideal atmosphere will be employed\n\n";
            }
            rtOpts.gradN = [&](const Vector3D& pos, const double& r){
                return this->getRefrac(pos, r);
            };
        }
    }
    

    void FlashMap::operator()()
    {
        setPosition(lon, lat);
        setTurbulence(turb);
        if(multi > 0 && turb && turbAvailable){
            std::cout << "Performing multiple map calculation\n";
            multiMap();
        }
        else singleMap();
    }



    
    void FlashMap::singleMap()
    {
        double dphi_t = 2*nb::pi/NTHREAD;
        std::vector<std::thread> threads(0);
        std::vector<fmap> map_t(NTHREAD);
        Vector2D rs;

        // calculate discretization of the atmosphere
        dr = dr_n*H*S/(20*N*R);
        dphi = dphi_n*2*nb::pi/(5*N);
        rtOpts.dt = dt_n*10*R/(R_std);

        // clear all auxiliary maps
        for(fmap& m : map_t){
            m.N = N;
            m.S = S;
            m.L = L;
            m.lat = lat;
            m.lon = lon;
            m.rtOpts = rtOpts;

            m.A = 0;
            m.n = 0;
            m.ntot = 0;
            m.clear();
            m.resize(N, std::vector<double>(N, 0));
        }

        // clear final map
        A = 0;
        n = 0;
        ntot = 0;
        clear();
        resize(N, std::vector<double>(N, 0));

        // print details
        printParameters();
        printFlashpars();

        std::cout << "Discretization:\n";
        std::cout << "dr: " << dr << " km\n";
        std::cout << "dphi: " << dphi << " rad\n";
        std::cout << "dt: " << rtOpts.dt << "\n\n";

        if(turb && turbAvailable){
            std::cout << "Turbulence is active.\n\n";
        }
        else{
            std::cout << "Turbulence is not active.\n\n";
        }

        auto start = std::chrono::high_resolution_clock::now();
        // calculate optimum position
        std::cout << "Calculating Optimum Position . . ." << std::endl;
        rs = getStartingRadius();
        std::cout << "Optimum at: " << rs << std::endl;
        if(rs[X] == R || rs[Y] == R){
            std::cout << "WARNING: Detector is totally or partially covered by the planet\n";
        }

        std::cout << "Starting Threads . . .\n" << std::endl;
        for(uint i = 0; i < NTHREAD; ++i){
            threads.push_back(std::thread([&, i]{this->mapThread<double>(map_t[i], i*dphi_t, (i+1)*dphi_t, rs, i+1); }));
        }

        for(std::thread& t : threads){
            t.join();
        }

        for(const fmap& m : map_t){
            (*this) += m;
        }

        normalize(A*sqr(N)/(n*sqr(S)));
        auto end = std::chrono::high_resolution_clock::now();
        finalData(end - start);
    }

    void FlashMap::multiMap()
    {
        std::string header = turbMap->getHeader(), filename;
        fmap multiAux;
        int k = 1;

        multiAux.N = N;
        multiAux.S = S;
        multiAux.L = L;
        multiAux.lat = lat;
        multiAux.lon = lon;
        multiAux.rtOpts = rtOpts;
        multiAux.A = 0;
        multiAux.n = 0;
        multiAux.ntot = 0;
        multiAux.clear();
        multiAux.resize(N, std::vector<double>(N, 0));

        for(; k <= multi; ++k){
            filename = header + '_' + std::to_string(k) + ".dat";
            if(!fs::exists(filename)) break;
            turbMap->readTurbFile(filename);
            std::cout << "In file " << filename << "...\n";
            singleMap();
            
            for(luint i = 0; i < size(); ++i){
                for(luint j = 0; j < size(); ++j){
                    multiAux[i][j] += at(i)[j] + at(j)[size() - 1 - i] + at(size() - 1 - i)[size() - 1 - j] + at(j)[size() - 1 - i];
                }
            }

            multiAux.A += A;
            multiAux.n += n;
            multiAux.ntot += ntot;
        }

        // clear final map
        A = 0;
        n = 0;
        ntot = 0;
        clear();
        resize(N, std::vector<double>(N, 0));
        (*this) += multiAux;
        normalize(1.0/(4.0*k));
    }


    template<typename U>
    FlashOutput FlashMap::propagateRay(const Vector3D& x0s, const Vector3D& v0s) const
    {
        Vector3D x0 = rtt*x0s, v0 = rtt*v0s;
        Vector2D xy;
        FlashOutput fOut;
        double phiB = 0;

        // cast ray to the atmosphere sphere
        if (std::is_same<U, complex>::value) {
            x0 = ray::ray2spherePhase(x0, v0, phiB, rtOpts.l, {R_atm, beta[Y], beta[Z]});
        }
        else {
            x0 = ray::ray2sphere(x0, v0, {R_atm, beta[Y], beta[Z]});
        }

        if(std::isnan(x0[X])){
            // cast ray to the detector
            xy = ray::rayZplane(x0s, v0s, L);

            // map to the detector plane
            map2detector(xy, fOut.i, fOut.j);
            
            fOut.T = 1;
            fOut.phi = 0;

            return fOut;
        }

        // propagate ray through the atmosphere
        auto [x, v, rmin, T, phi] = ray::rayTrace(x0, v0, rtOpts);
        fOut.T = T;

        // check if surface is hit
        if(T == 0) return fOut;

        // cast ray to the detector
        xy = ray::rayZplane(rti*x, rti*v, L);

        // map to the detector plane
        map2detector(xy, fOut.i, fOut.j);

        // get common phase if wave optics
        if  (std::is_same<U, complex>::value) {
            ray::ray2spherePhase(x, v, phi, rtOpts.l, {0.9*L,1,1},{xy[X], xy[Y], L});
            fOut.phi = phi;
        }
        
        return fOut;
    }

    /*FlashOutput FlashMap::propagateRayPhase(const Vector3D& x0s, const Vector3D& v0s) const
    {
        Vector3D x0 = rtt*x0s, v0 = rtt*v0s;
        Vector2D xy;
        FlashOutput fOut;
        double phiB = 0;

        x0 = ray::ray2spherePhase(x0, v0, phiB, rtOpts.l, {R_atm, beta[Y], beta[Z]});

        if(std::isnan(x0[X])){
            // cast ray to the detector
            xy = ray::rayZplane(x0s, v0s, L);

            // map to the detector plane
            map2detector(xy, fOut.i, fOut.j);
            
            fOut.T = 1;
            fOut.phi = 0;

            return fOut;
        }

        // propagate ray through the atmosphere
        auto [x, v, rmin, T, phi] = ray::rayTracePhase(x0, v0, phiB, rtOpts);
        fOut.T = T;

        // check if surface is hit
        if(T == 0) return fOut;

        // cast ray to the detector
        xy = ray::rayZplane(rti*x, rti*v, L);

        // map to the detector plane
        map2detector(xy, fOut.i, fOut.j);

        // get common phase
        ray::ray2spherePhase(x, v, phi, rtOpts.l, {0.9*L,1,1},{xy[X], xy[Y], L});
        fOut.phi = phi;

        return fOut;
    }*/




    void FlashMap::map2detector(const Vector2D& xy, int& i, int& j) const
    {
        j = std::floor(N*(xy[Y]/S + 0.5));
        i = std::floor(N*(0.5 - xy[X]/S));
    }



    Vector2D FlashMap::getStartingRadius()
    {
        Vector3D v0 = Vector3D({0,0,1});
        Vector2D p_s = {0, nb::pi/2}, res = {R,R}, r_s;
        double drn, R = this->R, b = beta[Y], t2t;
        bool tf, otherside = false;
        int aux_ij;

        if(lat != 0){
            // calculate properties of the cross-sectioned ellipse
            t2t = (lat == 90) ? 1/beta[Z] : std::sqrt((1.0 + sqr(std::tan(latP)))/(1.0 + sqr(beta[Z]*std::tan(latP))));
            R = this->R*t2t;
            b = beta[Y]*t2t;
        }

        if(turb && turbAvailable){
            rtOpts.gradN = [&](const Vector3D& pos, const double& r){
                return this->getRefrac(pos, r);
            };
        }

        // approximation techniques to obtain the initial height
        r_s[X] = 0.5*H*(1 - 0.5*H/R)*std::log(2*nb::pi*sqr(b*n0*L)/(R*H)) - 9*sqr(H)/(8*R*b);
        r_s[Y] = 0.5*H*(1 - 0.5*H/this->R)*std::log(2*nb::pi*sqr(beta[Y]*n0*L)/(this->R*H)) - 9*sqr(H)/(8*this->R*beta[Y]);

        // transformation to initial radial position
        r_s[X] = r_s[X]*(1 + 2*H/R) + R/b;
        r_s[Y] = r_s[Y]*(1 + 2*H/this->R) + this->R/beta[Y];

        for(luint k = 0; k < 2; ++k){
            Vector3D x0 = Vector3D({r_s[k]*std::cos(p_s[k]), r_s[k]*std::sin(p_s[k]), -2*R});
            auto [i0, j0, TT, PHI] = propagateRay<double>(x0, v0);

            // check location of the initial guess
            // in the detector plane, move to the outside
            if(i0 >= 0 && i0 < (int)N && j0 >= 0 && j0 < (int)N){
                drn = dr;
                tf = true;
            }
            // outside the detector plane, but on the opposite of the expected side
            else if((i0 >= (int)N && j0 >= 0) || (j0 < 0 && i0 >= 0)){
                drn = dr;
                tf = false;
                otherside = true;
            }
            // outisde the detector on the expected side
            else{
                drn = -dr;
                tf = false;
            }

            for(double r = r_s[k] + drn; r > R; r += drn){
                x0 = Vector3D({r*std::cos(p_s[k]), r*std::sin(p_s[k]), -2*R});
                auto [i,j,TT1,PHI1] = propagateRay<double>(x0, v0);

                if((i < 0 || j < 0 || i >= (int)N || j >= (int)N) == tf){
                    if(otherside){
                        otherside = false;
                        tf = !tf;
                        continue;
                    }
                    
                    res[k] = r;
                    break;
                }
                // dirty tricks to accelerate convergence
                else if((aux_ij = std::abs(i - (int)N/2) + std::abs(j - (int)N/2)) > (int)(1000*N)){
                    drn = 10000*((drn > 0) ? dr : -dr);
                }
                else if(aux_ij > (int)(100*N)){
                    drn = 1000*((drn > 0) ? dr : -dr);
                }
                else if(aux_ij > (int)N){
                    drn = 100*((drn > 0) ? dr : -dr);
                }
                else{
                    drn = (drn > 0) ? dr : - dr;
                }
            }
        }

        if(turb && turbAvailable){
            rtOpts.gradN = [&](const Vector3D& pos, const double& r){
                return this->getRefracTurb(pos, r);
            };

            //res[X] += 5*N*dr;
            //res[Y] += 5*N*dr;
        }

        return res;
    }

    inline int FlashMap::squareIntersect(const int& i1, const int& j1, const int& i2, const int& j2) const
    {
        Vector2D v = Vector2D({double(i2-i1), double(j2-j1)});
        double t, aux;

        if(v[X] != 0){
            t = ((int)N - i1)/v[X];
            aux = t*v[Y] + j1;
            if(aux >= 0 && aux <= N){
                return t <= 0;
            }

            t = -i1/v[X];
            aux = t*v[Y] + j1;
            if(aux >= 0 && aux <= N){
                return t <= 0;
            }
        }
        if(v[Y] != 0){
            t = ((int)N - j1)/v[Y];
            aux = t*v[X] + i1;
            if(aux >= 0 && aux <= N){
                return t <= 0;
            }

            t = -j1/v[Y];
            aux = t*v[X] + i1;
            if(aux >= 0 && aux <= N){
                return t <= 0;
            }
        }

        return -1;
    }

    template<typename U>
    inline void FlashMap::singleLoop(fmap<U>& m, const double& cp, const double& sp, const double& rs, const bool& down) const
    {
        double drn = (down) ? -dr : dr;
        double dro = drn;
        bool inside = false;
        int aux_ij;
        Vector3D x0;

        for (double r = rs; r > R; r += drn) {
            m.ntot += 1;
            x0 = Vector3D({r*cp, r*sp, -R_atm});
            auto [i, j, T, phi_ph] = propagateRay<U>(x0, v0);
            if(T == 0) break;
            else if (i < 0 || j < 0 || i >= (int)N || j >= (int)N) {
                if(inside) {
                    break;
                }
                continue;
                // dirty tricks to accelerate convergence
                if ((aux_ij = std::abs(i - (int)N/2) + std::abs(j - (int)N/2)) > (int)(1000*N)) {
                    drn = 10000*(dro);
                }
                else if (aux_ij > (int)(100*N)) {
                    drn = 1000*(dro);
                }
                else if (aux_ij > (int)N) {
                    drn = 100*(dro);
                }
                else{
                    drn = dro;
                }
            }
            inside = true;
            add2map(m, T, phi_ph, i, j);
            m.A += dr*r*dphi;
            m.n += 1;
        }
    }

    template<>
    inline void FlashMap::add2map<double>(fmap<double>& m, const double& T, const double& phi, const int& i, const int& j)
    {
        m[i][j] += T;
    }

    template<>
    inline void FlashMap::add2map<complex>(fmap<complex>& m, const double& T, const double& phi, const int& i, const int& j)
    {
        m[i][j] += T*std::exp(complex(0, phi));
    }

    template <typename U>
    void FlashMap::mapThread(fmap<U>& m, const double& phi_s, const double& phi_e, const Vector2D& rs, const int& iThread) const
    {
        double cp, sp, rss;
        uint icur = 0;
        int dir;
        //Vector3D x0;
        //bool inside;
        if(m.size() != N){
            m.clear();
            m.resize(N, std::vector<double>(N, 0));
        }

        for(double phi = phi_s; phi < phi_e; phi += dphi){
            cp = std::cos(phi);
            sp = std::sin(phi);
            rss = (rs[X] - rs[Y])*sqr(cp) + rs[Y];
            if(100*(phi - phi_s)/(phi_e - phi_s) >= icur){
                printf("Thread %d: at %d%%\n", iThread, icur);
                icur += 10;
            }
            m.ntot += 1;
            auto [i, j, T, phi_p] = propagateRay<double>(Vector3D({rss*cp, rss*sp, -R_atm}), v0);
            if(i < 0 || j < 0 || i >= (int)N || j >= (int)N){
                //LOG
                m.ntot += 1;
                auto [i2, j2, T, phi_p] = propagateRay<double>(Vector3D({(rss+30*dr)*cp, (rss + 30*dr)*sp, -R_atm}), v0);

                if(i2 < 0 || j2 < 0 || i2 >= (int)N || j2 >= (int)N){
                    if((dir = squareIntersect(i, j, i2, j2)) == -1) continue;
                    singleLoop<U>(m, cp, sp, rss + ((dir) ? -dr : dr), dir);
                }
                else{
                    singleLoop<U>(m, cp, sp, rss+dr, false);
                }
            }
            else{
                singleLoop<U>(m, cp, sp, rss-dr, true);
                singleLoop<U>(m, cp, sp, rss, false);
            }
            
        }
    }

    template class fmap<double>;
    template class fmap<complex>;
};