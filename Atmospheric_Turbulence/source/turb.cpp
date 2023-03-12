#include "turb.hpp"

Turbulence::Turbulence(const turbOpts& topts)
{
    Co2 = (4*std::sqrt(M_PI)*G16*std::pow(topts.M, 4)/(45*G53))*std::pow(2*M_PI*topts.R/topts.l0, 2.0/3.0);

    std::cout << "Co2 = " << Co2 << std::endl;
    std::cout << "Co2_old = " << (12.0*M_PI*std::pow(topts.M, 4)/5.0)*std::pow(2*M_PI*topts.R/topts.l0, 2.0/3.0) << std::endl;
}

void Turbulence::genCoef(const int& lplus, const int& N)
{
    double aux_e;
    std::random_device dev;
    std::mt19937 gen(dev());
    std::normal_distribution<double> dist;

    lmax_e = lplus;

    Clm.clear();
    Clm.resize(lmax_e + 1);
    E.clear();
    Er.clear();

    std::cout << "Generating spherical harmonics . . ." << std::endl;
    if(lmax_e > lmax) calcAll(lmax_e, N);
    std::cout << "Done." << std::endl;
    

    std::cout << "Generating refractivity coefficients . . ." << std::endl;
    for(int ll = 1; ll <= lmax_e; ++ll){
        aux_e = 0;
        Clm[ll].clear();
        Clm[ll].resize(ll+1);

        dist = std::normal_distribution<double>(0, std::sqrt(Co2/std::pow(ll, 8.0/3.0)));
        Er.emplace_back(ll, Co2*std::pow(ll, 8.0/3.0));
        
        for(int m = 0; m <= ll; ++m){
            Clm[ll][m] = {dist(gen), dist(gen)};
            aux_e += Clm[ll][m][0]*Clm[ll][m][0] + Clm[ll][m][1]*Clm[ll][m][1];
        }
        E.emplace_back(ll, aux_e);
    }
    std::cout << "Done." << std::endl;

    isGen = false;
}

void Turbulence::genMap(const int& Nphi, const double& tet_max)
{
    double dp = 2*M_PI/(Nphi-1);//, phi, tet;
    int Ntet = std::floor((tet_max/M_PI)*Nphi), is = 0, ie;
    std::vector<std::thread> th;

    th.clear();
    map.clear();
    g_tet.clear();
    g_phi.clear();

    map.resize(Nphi);
    g_tet.resize(Nphi);
    g_phi.resize(Nphi);

    std::cout << "Setting up threads . . ." << std::endl;
    for(uint i = 0; i < nthread; ++i){
        if(i != nthread - 1) ie = (i+1)*Nphi/nthread;
        else ie = Nphi;
        th.push_back(std::thread([this, is, ie, Ntet, dp] { this->genThread(is, ie, Ntet, dp); }));
        is = ie;
    }

    std::cout << "Generating atmosphere . . ." << std::endl;
    for(uint i = 0; i < nthread; ++i){
        th[i].join();
        std::cout << "Joined thread " << i+1 << std::endl;
    }
    std::cout << "Done." << std::endl;

    /*for(int i = 0; i < Nphi; ++i){
        map[i].resize(Ntet, 1);
        g_tet[i].resize(Ntet, 0);
        g_phi[i].resize(Ntet, 0);
        phi = dp*i;
        if(i % 10 == 0) std::cout << i << std::endl;
        for(int j = 0; j < Ntet; ++j){
            tet = dp*(j - (Ntet-1)/2.0) + M_PI_2;
            //if(i == j) std::cout << tet << std::endl;
            for(int l = 0; l <= lmax_e; ++l){
                map[i][j] += Clm[l][0][0]*(*this)(l, 0, tet, phi);
                g_tet[i][j] += Clm[l][0][0]*gtet(l, 0, tet, phi);
                g_phi[i][j] += Clm[l][0][0]*gphi(l, 0, tet, phi);
                for(int m = 1; m <= l; ++m){
                    map[i][j] += Clm[l][m][0]*(*this)(l, m, tet, phi) + Clm[l][m][1]*(*this)(l, -m, tet, phi);
                    g_tet[i][j] += Clm[l][m][0]*gtet(l, m, tet, phi) + Clm[l][m][1]*gtet(l, -m, tet, phi);
                    g_phi[i][j] += Clm[l][m][0]*gphi(l, m, tet, phi) + Clm[l][m][1]*gphi(l, -m, tet, phi);
                }
            }
        }
    }
    //max_tet = tet - M_PI_2;*/
    max_tet = dp*(Ntet-1)/2.0;
    isGen = true;
}

void Turbulence::genThread(const int& is, const int& ie, const int& Ntet, const double& dp)
{
    double tet, phi;

    for(int i = is; i < ie; ++i){
        map[i].resize(Ntet, 1);
        g_tet[i].resize(Ntet, 0);
        g_phi[i].resize(Ntet, 0);
        phi = dp*i;
        //if(i % 10 == 0) std::cout << i << std::endl;
        for(int j = 0; j < Ntet; ++j){
            tet = dp*(j - (Ntet-1)/2.0) + M_PI_2;
            //if(i == j) std::cout << tet << std::endl;
            for(int l = 1; l <= lmax_e; ++l){
                map[i][j] += Clm[l][0][0]*(*this)(l, 0, tet, phi);
                g_tet[i][j] += Clm[l][0][0]*gtet(l, 0, tet, phi);
                g_phi[i][j] += Clm[l][0][0]*gphi(l, 0, tet, phi);
                for(int m = 1; m <= l; ++m){
                    map[i][j] += Clm[l][m][0]*(*this)(l, m, tet, phi) + Clm[l][m][1]*(*this)(l, -m, tet, phi);
                    g_tet[i][j] += Clm[l][m][0]*gtet(l, m, tet, phi) + Clm[l][m][1]*gtet(l, -m, tet, phi);
                    g_phi[i][j] += Clm[l][m][0]*gphi(l, m, tet, phi) + Clm[l][m][1]*gphi(l, -m, tet, phi);
                }
            }
        }
    }
}

void Turbulence::printMap(const std::string& filename)
{
    std::fstream fp;

    if(!isGen){
        std::cout << "No map has been generated with current coefficients\n" << std::endl;
        return;
    }

    fp.open(filename, std::ofstream::out);
    fp << max_tet << " " << map.size() << " " << map[0].size() << std::endl;
    for(uint i = 0; i < map.size(); ++i){
        for(uint j = 0; j < map[i].size(); ++j){
            fp << map[i][j] << " " << g_tet[i][j] << " " << g_phi[i][j] << std::endl;
        }
    }
    fp.close();
}

void Turbulence::plotMap(const std::string& filename, const int& div)
{
    Gnuplot gp;
    std::ostringstream ss1, ss2;
    int nmax = std::floor(max_tet*div/M_PI);

    if(isGen){
        std::fstream fp("aux.dat", std::ofstream::out);
        for(uint i = 0; i < map[0].size(); ++i){
            for(uint j = 0; j < map.size(); ++j){
                fp << 100.0*(map[j][i] - 1) << " ";
            }
            fp << std::endl;
        }
        fp.close();
    }
    else return;
    
    

    ss1 << "($1*(2*pi/ " << map.size() << "))";
    ss2 << "($2*(" << 2*max_tet/map[0].size() << ") - " << max_tet <<")";
    
    
    gp << "reset\nunset key\n";
    if(filename == ""){
        gp << "set terminal wxt size 800,600 enhanced font 'Verdana,8' persist\n";
    
    }
    else{
        gp << "set terminal png size 1600,400 enhanced font 'Verdana,15' enhanced\nset output '" + filename + "'\n";
    }
    gp << "set encoding utf8\nset termoption enhanced\n";
    gp << "set style line 11 lc rgb '#808080' lt 1\n";
    gp << "set border 3 front ls 11\n";
    gp << "set tics nomirror out scale 0.75\n";
    gp << "set xtics pi\nset ytics pi\n";
    gp << "set xtics (\"0\" 0, \"π\" pi, \"2π\" 2*pi)\n";
    gp << "set palette defined (0 0 0 0, 1 1 1 1)\n";
    gp << "set yrange [" + std::to_string(-max_tet) + ":" + std::to_string(max_tet) + "]\nset ylabel 'θ - rad' enhanced\n";
    gp << "set ytics (";
    for(int n = -nmax; n <= nmax; ++n){
        if(n != -nmax) gp << ", ";

        if(std::abs(n) > 1){
            gp << "\"" + std::to_string(n) + "π/" + std::to_string(div);
        }
        else{
            if(n == 0) gp << "\"0";
            else{
                if(n == 1) gp << "\"π/" + std::to_string(div);
                else gp << "\"-π/" + std::to_string(div);
            } 
        }
        gp << "\" " + std::to_string(n) + "*pi/" + std::to_string(div);
    }
    gp << ")\n";
    gp << "set xrange [0:2*pi]\nset xlabel 'φ - rad' enhanced\n";
    gp << "set view map\nset size ratio -1\n";
    gp << "set pm3d interpolate 2,2\nsplot 'aux.dat' matrix using " + ss1.str() + ":" + ss2.str() + ":($3) with pm3d\n";
}