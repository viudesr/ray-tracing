#ifndef VECTOR_HPP
#define VECTOR_HPP

#include <random>
#include <math.h>
#include <omp.h>

static inline double sqr(double x) { return x * x; }

std::default_random_engine engine[12];
std::uniform_real_distribution<double> uniform;

class Vector {
public:
    Vector(double x = 0, double y = 0, double z = 0) {
        x = x;
        y = y;
        z = z;
    }
    double norm2() {
        return sqr(x) + sqr(y) + sqr(z);
    }
    void normalize() {
        double norm = sqrt(this->norm2());
        x /= norm;
        y /= norm;
        z /= norm;
    }
    double x;
    double y;
    double z;
};
double dot(Vector& A, Vector& B) {
    return A.x * B.x + A.y * B.y + A.z * B.z;
};
#endif