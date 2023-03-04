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
        this->x = x;
        this->y = y;
        this->z = z;
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

Vector operator+(const Vector & A, Vector& B) {
    return Vector(A.x + B.x, A.y + B.y, A.z + B.z);
}
Vector operator-(const Vector & A, Vector& B) {
    return Vector(A.x - B.x, A.y - B.y, A.z - B.z);
}
Vector operator*(const Vector & A, Vector& B) {
    return Vector(A.x * B.x, A.y * B.y, A.z * B.z);
}
double dot(Vector& A, Vector& B) {
    return A.x * B.x + A.y * B.y + A.z * B.z;
};
#endif