#ifndef UTILS_HPP
#define UTILS_HPP

#define _USE_MATH_DEFINES

#include <cmath>
#include <complex>
#include <array>
#include <vector>
#include <numeric>
#include <numbers>
#include <ranges>
#include <functional>
#include <algorithm>
#include <string>
#include <limits>
#include <tuple>

#include <iostream>
#include <fstream>
#include <filesystem>

// Useful for debugging
#define LOG {std::cout << "IN LINE " << __LINE__ << " OF FILE " << __FILE__ << std::endl;}
#ifndef NAN
    #define NAN std::nan("")
#endif

// namespaces
namespace rv = std::views;
namespace fs = std::filesystem;
namespace nb = std::numbers;

namespace std::numbers {
    constexpr double inf = std::numeric_limits<double>::infinity();
};

// Positional defines
#define X 0
#define Y 1
#define Z 2

#define BX2 0
#define BZ2 1
#define BXZ 2

constexpr double rad2deg = 180*nb::inv_pi;
constexpr double deg2rad = nb::pi/180;

// Useful datatypes
using luint = long unsigned int;
using complex = std::complex<double>;

template<typename T, luint N>
using Vector = std::array<T,N>;
using Vector3D = std::array<double, 3>;
using Vector2D = std::array<double, 2>;

template<typename T, luint N, luint M>
using Matrix = std::array<std::array<T,M>,N>;
using Matrix3D = std::array<std::array<double, 3>, 3>;

using pairXY = std::pair<double, double>;
using VpairXY = std::vector<pairXY>;
template<typename T>
using Vmap = std::vector<std::vector<T>>;

// useful function definitions
using funcScaSca = std::function<double(const double&)>;
using funcScaVec = std::function<double(const Vector3D&)>;
using funcScaVecSca = std::function<double(const Vector3D&, const double&)>;
using funcVecVecSca = std::function<Vector3D(const Vector3D&, const double&)>;
using funcExit = std::function<bool(const Vector3D&, const double&)>;

// Specific refractivity function and output
struct refracOutput {
    double n1;
    Vector3D gradN;
};
using funcRefrac = std::function<refracOutput(const Vector3D&, const double&)>;

constexpr double t_95[] = {1.96, 12.706, 4.303, 3.182, 2.776, 2.571, 2.447, 2.365, 2.306, 2.262, 2.228, 2.201, 2.179, 2.160, 2.145, 2.131, 2.120, 2.110, 2.101, 2.093, 2.086, 2.080, 2.074, 2.069, 2.064, 2.060, 2.056, 2.052, 2.048, 2.045, 2.042};

/*
    Definition of vector operators
    - Output
    - Addition
    - Subtraction
    - Multiplication
    - Division
    - Inner and outer product
*/

template<typename T, luint N>
inline std::ostream& operator<<(std::ostream& os, const Vector<T,N>& v)
{
    os << "[ ";
    for(const T& elem : v){
        os << elem << ' ';
    }
    os << ']';

    return os;
}

template<typename T, luint N>
inline T operator*(const Vector<T,N>& l, const Vector<T,N>& r)
{
    T res = 0;

    for(luint i = 0; i < N; ++i){
        res += l[i]*r[i];
    }

    return res;
}

template<typename T, luint N>
inline Vector<T,N> operator+(const Vector<T,N>& l, const Vector<T,N>& r)
{
    Vector<T,N> res;

    for(luint i = 0; i < N; ++i){
        res[i] = l[i] + r[i];
    }

    return res;
}

template<typename T, luint N>
inline Vector<T,N> operator+=(Vector<T,N>& l, const Vector<T,N>& r)
{
    for(luint i = 0; i < N; ++i){
        l[i] += r[i];
    }

    return l;
}

template<typename T, luint N>
inline Vector<T,N> operator-(const Vector<T,N>& l, const Vector<T,N>& r)
{
    Vector<T,N> res;

    for(luint i = 0; i < N; ++i){
        res[i] = l[i] - r[i];
    }

    return res;
}

template<typename T, luint N>
inline Vector<T,N> operator-=(const Vector<T,N>& l, const Vector<T,N>& r)
{
    for(luint i = 0; i < N; ++i){
        l[i] -= r[i];
    }

    return l;
}

template<typename T, luint N>
inline Vector<T,N> operator*(const Vector<T,N>& l, const T& r)
{
    Vector<T,N> res;

    for(luint i = 0; i < N; ++i){
        res[i] = l[i]*r;
    }

    return res;
}

template<typename T, luint N>
inline Vector<T,N> operator*(const T& l, const Vector<T,N>& r)
{
    return r*l;
}

template<typename T, luint N>
inline Vector<T,N> operator*=(Vector<T,N>& l, const T& r)
{
    for(T& elem : l){
        elem *= r;
    }

    return l;
}

template<typename T, luint N>
inline Vector<T,N> operator/(const Vector<T,N>& l, const T& r)
{
    Vector<T,N> res;

    for(luint i = 0; i < N; ++i){
        res[i] = l[i]/r;
    }

    return res;
}

template<typename T, luint N>
inline Vector<T,N> operator/=(Vector<T,N>& l, const T& r)
{
    for(T& elem : l){
        elem /= r;
    }

    return l;
}

template<typename T, luint N>
inline Vector<T,N> operator^(const Vector<T,N>& l, const Vector<T,N>& r)
{
    Vector<T,N> res;

    for(luint i = 0; i < N; ++i){
        res[i] = l[i]*r[i];
    }

    return res;
}

template<typename T, luint N>
inline double normP2(const Vector<T,N>& v)
{
    return std::sqrt(v*v);
}

template<typename T, luint N>
inline double normPn(const Vector<T,N>& v, const double& p)
{
    double res = 0;

    for(const T& elem : v){
        res += std::pow(elem, p);
    }

    return std::pow(res, 1/p);
}

inline Vector3D cross(const Vector3D& l, const Vector3D& r)
{
    Vector3D res;

    res[X] = l[Y]*r[Z] - l[Z]*r[Y];
    res[Y] = l[Z]*r[X] - l[X]*r[Z];
    res[Z] = l[X]*r[Y] - l[Y]*r[Z];

    return res;
}



/*
    Definition of matrix operators
    - Output
    - Matrix multiplication
    - Vector-Matrix multiplication
*/

template<typename T, luint N, luint M>
inline std::ostream& operator<<(std::ostream& os, Matrix<T,N,M>& m)
{
    os << "[\n";
    for(const Vector<T,M>& v : m){
        os << v << '\n';
    }
    os << ']';

    return os;
}

template<typename T, luint N, luint K, luint M>
inline Matrix<T,N,M> operator*(const Matrix<T,N,K>& l, const Matrix<T,K,M>& r)
{
    Matrix<T,N,M>& res;

    for(luint i = 0; i < N; ++i){
        for(luint j = 0; j < M; ++j){
            res[i][j] = 0;
            for(luint k = 0; k < K; ++k){
                res[i][j] += l[i][k]*r[k][j];
            }
        }
    }

    return res;
}

template<typename T, luint N, luint M>
inline Vector<T,N> operator*(const Matrix<T,N,M>& l, const Vector<T,M>& r)
{
    Vector<T,N> res;

    for(luint i = 0; i < N; ++i){
        res[i] = l[i]*r;
    }

    return res;
}

template<typename T, luint N, luint M>
inline Vector<T,M> operator*(const Vector<T,N>& l, const Matrix<T,N,M>& r)
{
    Vector<T,M> res;

    for(luint j = 0; j < M; ++j){
        res[j] = 0;
        for(luint i = 0; i < N; ++i){
            res[j] += l[i]*r[i][j];
        }
    }

    return res;
}



// resolving conflicts
// defition required due to conflicting templates
inline Matrix3D operator*(const Matrix3D& l, const Matrix3D& r)
{
    Matrix3D res;

    for(luint i = 0; i < 3; ++i){
        for(luint j = 0; j < 3; ++j){
            res[i][j] = 0;
            for(luint k = 0; k < 3; ++k){
                res[i][j] += l[i][k]*r[k][j];
            }
        }
    }

    return res;
}

inline Vector3D operator*(const Vector3D& l, const Matrix3D& r)
{
    Vector3D res;

    for(luint j = 0; j < 3; ++j){
        res[j] = 0;
        for(luint i = 0; i < 3; ++i){
            res[j] += l[i]*r[i][j];
        }
    }

    return res;
}


/*
    Useful functions
    - Squaring
    - Average
    - Variance
*/

template<typename T>
inline T sqr(const T& x)
{
    return x*x;
}

template<typename T>
inline T avg(const std::vector<T>& v)
{
    return (v.size()) ? std::accumulate(v.begin(), v.end(), T(0))/v.size() : T(0);
}

template<typename T>
inline T var(const std::vector<T>& v)
{
    if(v.size() == 0) return T(0);

    std::vector<T> aux;
    T aux_avg = avg(v);

    for(const T& elem : v){
        aux.emplace_back(sqr(elem - aux_avg));
    }

    return avg(aux);
}

template<typename T>
inline T stdDev(const std::vector<T>& v)
{
    return std::sqrt(var(v));
}

inline double tStudentCritical95(const uint& n)
{
    /*if(n > 30){
        return t_95[0];
    }*/

    return t_95[(n > 30) ? 0 : n];
}

template<typename T> inline constexpr
int signum(T x, std::false_type is_signed)
{
    return T(0) < x;
}

template<typename T> inline constexpr
int signum(T x, std::true_type is_signed)
{
    return (T(0) < x) - (T(0) > x);
}

template<typename T> inline constexpr
int signum(T x)
{
    return signum(x, std::is_signed<T>());
}


#endif