#ifndef VECTOR_HPP
#define VECTOR_HPP

#include <random>
#include <math.h>
#include <omp.h>
#include <vector>

static inline double sqr(double x) { return x * x; }

std::default_random_engine engine[12]; // 12 threads
std::uniform_real_distribution<double> uniform;

class Vector {
public:
    Vector(double x = 0, double y = 0, double z = 0) {
        this->x = x;
        this->y = y;
        this->z = z;
    }

    const double& operator [] (int i) const {
        if (i==0) {
            return this->x;
        }
        if (i==1) {
            return this->y;
        }
        if (i==2) {
            return this->z;
        }
    }
    Vector& operator += (const Vector& A) {
        x += A[0];
        y += A[1];
        z += A[2];
        return *this;
    }
    Vector& operator += (const double& a) {
        x += a;
        y += a;
        z += a;
        return *this;
    }
    Vector& operator -= (const Vector& A) {
        x -= A.x;
        y -= A.y;
        z -= A.z;
        return *this;
    }
    Vector& operator -= (const double& a) {
        x -= a;
        y -= a;
        z -= a;
        return *this;
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

private:
    double x;
    double y;
    double z;
};

Vector operator+(const Vector &A, const Vector& B) {
    return Vector(A[0] + B[0], A[1] + B[1], A[2] + B[2]);
}
Vector operator+(const Vector &A, double& b) {
    return Vector(A[0] + b, A[1] + b, A[2] + b);
}
Vector operator+(double& b, const Vector &A) {
    return Vector(A[0] + b, A[1] + b, A[2] + b);
}
Vector operator-(const Vector &A, const Vector& B) {
    return Vector(A[0] - B[0], A[1] - B[1], A[2] - B[2]);
}
Vector operator-(const Vector &A, double& b) {
    return Vector(A[0] - b, A[1] - b, A[2] - b);
}
Vector operator-(double& b, const Vector &A) {
    return Vector(b - A[0], b - A[1], b - A[2]);
}
Vector operator-(const Vector& A) {
        return Vector(-A[0], -A[1], -A[2]);
    }
Vector operator*(const Vector &A, const Vector& B) {
    return Vector(A[0] * B[0], A[1] * B[1], A[2] * B[2]);
}
Vector operator*(const Vector &A, double& b) {
    return Vector(A[0] * b, A[1] * b, A[2] * b);
}
Vector operator*(double& b, const Vector &A) {
    return Vector(A[0] * b, A[1] * b, A[2] * b);
}
Vector operator/(const Vector &A, double b) {
    return Vector(A[0] / b, A[1] / b, A[2] / b);
}
Vector operator/(double &b, const Vector &A) {
    return Vector(b / A[0], b / A[1], b / A[2]);
}
double dot(const Vector& A, const Vector& B) {
    return A[0] * B[0] + A[1] * B[1] + A[2] * B[2];
};
Vector cross(const Vector& A, const Vector& B) {
    return Vector(A[1] * B[2] - B[1] * A[2], A[2] * B[0] - B[2] * A[0], A[0] * B[1] - B[0] * A[1]);
}
#endif