#include "test.hpp"

void plotLegPoly(const int& l, const int& m)
{
    PolyCalc p;
    std::vector<VpairXY> res(l-m+1);
    Gnuplot gp;
    auto plt = gp.plotGroup();
    int pol;

    p.calcAll(l);

    for(double x = -1; x <= 1.001; x += 0.01){  
        p.calcPoly(l, m, x);
        for(int i = m; i <= l; ++i){
            //std::cout << p.res[m][i] << std::endl;;
            res.at(i-m).emplace_back(x, p.normCoef[l][i]*p.calcPoly(l, i, x));
            //std::cout << p.res[m][i] << " ";
        }
    }

    pol = m;
    for(VpairXY& res1 : res){
        plt.add_plot1d(res1, "with lines title 'm = " + std::to_string(pol++) + "'");
    }

    gp << plt;
}


void plotLs(const int& l, const int& ngraph)
{
    VpairXY res;
    PolyCalc p;
    p.calcAll(l);

    for(int m = 0; m <= l; ++m){
        Gnuplot gp;
        auto plt = gp.plotGroup();
        for(int ll = m; (ll <= l) && (ll-m <= ngraph); ++ll){
            res = p.poly.at(ll).at(m).getVals();
            plt.add_plot1d(res, "with lines title 'l = " + std::to_string(ll) + "'");
        }
        gp << plt;
    }
}

void plotSphereHarm(const int& l, const int& m, const int& N)
{
    Gnuplot gp;
    auto plt = gp.splotGroup();
    std::vector<std::vector<double>> x(N), y(N), z(N);
    SphereHarm p;
    double dtet = M_PI/(N+1), dphi = 2*M_PI/(N-1);
    double tet, phi, r;

    p.calcAll(l);

    for(int i = 0; i < N; ++i){
        tet = dtet*(i+1);
        x[i].resize(N);
        y[i].resize(N);
        z[i].resize(N);
        for(int j = 0; j < N; ++j){
            phi = dphi*j;
            r = std::abs(p(l, m, tet, phi));
            x[i][j] = r*std::sin(tet)*std::cos(phi);
            y[i][j] = r*std::sin(tet)*std::sin(phi);
            z[i][j] = r*std::cos(tet);
        }
    }

    //gp << "set dgrid3d 50,50,1\n";
    //gp << "set hidden3d\n";

    gp << "set xrange [-1:1]\n";
    gp << "set yrange [-1:1]\n";
    gp << "set zrange [-1:1]\n";
    plt.add_plot2d(std::make_tuple(x, y, z), "with lines");
    gp << plt;
    getchar();
}

void plotCn2()
{
    Gnuplot gp;
    auto plot = gp.plotGroup();
    VpairXY xy1, xy2;
    std::function<double(double)> Cn2_f;
    double H = 8, Cn2_0 = 2.0*std::pow(0.000275*(0.0309391487*0.0309391487),2)/(5.0 * std::pow(100000/(2*M_PI), 2.0/3.0));

    Cn2_f = [](double h){
        double aux;
        if(h < 2.13) aux = -10.7025-4.3507*h+0.8141*h*h;
        else if(h < 10.34) aux = -16.2897+0.0335*h-0.0134*h*h;
        else aux = -17.0577-0.0449*h-0.0005*h*h+0.6181*exp(-0.5*std::pow((h-15.5617)/3.4666, 2.0));

        return std::pow(10.0, aux);
    };

    for(double h = 0; h <= 100; h += 0.1){
        xy1.emplace_back(Cn2_f(h), h);
        xy2.emplace_back(Cn2_0*exp(-2*h/H), h);
    }
    plot.add_plot1d(xy1, "with lines title 'Fit to data'");
    plot.add_plot1d(xy2, "with lines title 'Simple model'");

    gp << "set logscale x\n";
    gp << "set xlabel 'C_n^2 - m^{-2/3}'\n";
    gp << "set ylabel 'Altitude h - km'\n";
    gp << plot;
}

#if false

void plotLsSphere(const int& l, const int& ngraph)
{
    VpairXY res;
    SphereHarm p;

    for(int m = 0; m <= l; ++m){
        Gnuplot gp;
        auto plt = gp.plotGroup();
        for(int ll = m; (ll <= l) && (ll-m <= ngraph); ++ll){
            if(ll == 0 && m == 0) continue;
            res.clear();
            for(double x = -1; x <= 1.00001; x += 0.001){
                res.emplace_back(x, p.poly[ll][m](x));
            }
            plt.add_plot1d(res, "with lines title 'l = " + std::to_string(ll) + "'");
        }
        gp << plt;
    }
}

void plotMap(const Turbulence& t)
{
    Gnuplot gp;
    std::ostringstream ss1, ss2;

    if(t.isGen){
        std::fstream fp("aux.dat", std::ofstream::out);
        for(uint i = 0; i < t.map[0].size(); ++i){
            for(uint j = 0; j < t.map.size(); ++j){
                fp << 100.0*(t.map[j][i] - 1) << " ";
            }
            fp << std::endl;
        }
        fp.close();
    }
    
    

    ss1 << "($1*(2*pi/ " << t.map.size() << "))";
    ss2 << "($2*(" << 2*t.max_tet/t.map[0].size() << ") - " << t.max_tet <<")";
    
    

    gp << "reset\nunset key\n";
    gp << "set terminal wxt size 800,600 enhanced font 'Verdana,8' persist\n";
    gp << "set style line 11 lc rgb '#808080' lt 1\n";
    gp << "set border 3 front ls 11\n";
    gp << "set tics nomirror out scale 0.75\n";
    gp << "set palette defined (0 0 0 0, 1 1 1 1)\n";
    gp << "set yrange [" + std::to_string(-t.max_tet) + ":" + std::to_string(t.max_tet) + "]\nset ylabel '\\theta - rad'\n";
    gp << "set xrange [0:2*pi]\nset xlabel '\\phi - rad'\n";
    gp << "set view map\nset size ratio -1\n";
    gp << "set pm3d interpolate 2,2\nsplot 'aux.dat' matrix using " + ss1.str() + ":" + ss2.str() + ":($3) with pm3d\n";
}

inline double hav(const double& x)
{
    return 0.5*(1 - cos(x));
}

void monteCarloCn2(const Turbulence& t, const uint& n)
{
    if(!t.isGen) return;

    Gnuplot gp;
    auto plot = gp.plotGroup();
    std::mt19937 gen(std::random_device{}());
    std::uniform_real_distribution<double> dist(0,1), dist2(-7, -2);
    double tet1, tet2, phi1, phi2, h, n1, n2, minr = 100000, maxr = 0;
    double dtet = 2*t.max_tet/(t.map[0].size()-1), dphi = 2*M_PI/(t.map.size() - 1), dxy;
    VpairXY xy, xy_l;
    uint ii, jj;

    dxy = 1/(dtet*dphi);

    std::cout << 2*std::asin(std::sqrt(hav(M_PI))) << std::endl;

    for(uint i = 0; i < n; ++i){
        tet1 = t.max_tet*(2*dist(gen) - 1);
        phi1 = 2*M_PI*dist(gen);

        //tet2 = t.max_tet*(2*dist(gen) - 1);
        tet2 = tet1;
        phi2 = std::fmod(phi1 + std::pow(10, dist2(gen)), 2*M_PI);

        h = hav(tet2 - tet1) + cos(tet1)*cos(tet2)*hav(phi2 - phi1);
        h = 2*std::asin(std::sqrt(h));
        
        ii = std::floor(phi1/(2*M_PI)*(t.map.size()-1));
        jj = std::floor((tet1 + t.max_tet)/(2*t.max_tet)*(t.map[0].size()-1));
        if(ii == (t.map.size() - 1)) ii = 0;
        n1 = dxy*(t.map[ii][jj]*((jj+1)*dtet - tet1)*((ii+1)*dphi - phi1)
                + t.map[ii][jj+1]*(tet1 - jj*dtet)*((ii+1)*dphi - phi1)
                + t.map[ii+1][jj]*((jj+1)*dtet - tet1)*(phi1 - ii*dphi)
                + t.map[ii+1][jj+1]*(tet1 - jj*dtet)*(phi1 - ii*dphi));

        ii = std::floor(phi2/(2*M_PI)*(t.map.size()-1));
        jj = std::floor((tet2 + t.max_tet)/(2*t.max_tet)*(t.map[0].size()-1));
        if(ii == (t.map.size() - 1)) ii = 0;
        n2 = dxy*(t.map[ii][jj]*((jj+1)*dtet - tet2)*((ii+1)*dphi - phi2)
                + t.map[ii][jj+1]*(tet2 - jj*dtet)*((ii+1)*dphi - phi2)
                + t.map[ii+1][jj]*((jj+1)*dtet - tet2)*(phi2 - ii*dphi)
                + t.map[ii+1][jj+1]*(tet2 - jj*dtet)*(phi2 - ii*dphi));
    
        //if(std::isnan(h)) h = M_PI;
        xy.emplace_back(h, (n2-n1)*(n2-n1));
        xy_l.emplace_back(std::log(h), std::log((n2-n1)*(n2-n1)));
        if(h < minr) minr = h;
        else if(h > maxr) maxr = h;

        //std::cout << n1 << " " << n2 << std::endl;
        //std::cout << xy.back().second; getchar();
    }

    auto [a, b] = linReg(xy_l);
    std::cout << "Power law: " << a << std::endl;
    std::cout << "Cn2: " << exp(b) << ", Co2: " << ((t.isNew) ? t.Co2New : t.Co2) << std::endl;
    std::cout << maxr << " " << minr << std::endl;

    plot.add_plot1d(xy, "with dots");
    gp << "set logscale x\nset logscale y\n";
    gp << plot;
}






std::array<double,2> linReg(const VpairXY& xy)
{
    int n_p = 0;
    double a = 0, b = 0;

    for(uint i = 0; i < xy.size(); ++i){
        for(uint j = i+1; j < xy.size(); ++j){
            if(xy[i].first == xy[j].first) continue;
            ++n_p;
            a += (xy[i].second - xy[j].second)/(xy[i].first - xy[j].first);
        }
    }
    a /= n_p;

    for(uint i = 0; i < xy.size(); ++i){
        b += xy[i].second - a*xy[i].first;
    }
    b /= xy.size();

    return {a,b};
}

void testLinReg(const double& m, const double& b, const double& e)
{
    std::mt19937 gen(std::random_device{}());
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    VpairXY xy;

    for(double x = 0; x <= 20; x += 0.001)
    {
        xy.emplace_back(x, m*x + b + e*dist(gen));
    }

    auto [m_l, b_l] = linReg(xy);
    std::cout << m << " " << m_l << std::endl;
    std::cout << b << " " << b_l << std::endl;
}

#endif

void plotEnergy(const int& l, const int& n)
{

    Turbulence t;

    for(int i = 0; i < n; ++i){
        Gnuplot gp;
        auto plot = gp.plotGroup();
        t.genCoef(l);
        plot.add_plot1d(t.E, "with lines title 'Real fluctuations'");
        plot.add_plot1d(t.Er, "with lines title 'Theoretical fluctuations'");

        gp << "set terminal png size 1600,1200 enhanced font 'Verdana,30'\nset termoption enhanced\n";
        gp << "set output '../res/img_CN_" + std::to_string(i) + ".png'\nset key bottom right\n";
        gp << "set xlabel 'Length scale - l'\nset logscale x\n";
        gp << "set ylabel 'Fluctuation - n'\nset logscale y\n";
        gp << plot;
    }
    
}