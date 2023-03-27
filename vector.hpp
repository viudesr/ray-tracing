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
    explicit Vector(double x = 0, double y = 0, double z = 0) {
        coords[0] = x;
        coords[1] = y;
        coords[2] = z;
    }

    const double& operator [] (int i) const {
        return coords[i];
    }

    double& operator [] (int i) {
        return coords[i];
    }

    Vector& operator += (const Vector& A) {
        coords[0] += A[0];
        coords[1] += A[1];
        coords[2] += A[2];
        return *this;
    }

    Vector& operator += (const double& a) {
        coords[0] += a;
        coords[1] += a;
        coords[2] += a;
        return *this;
    }
    Vector& operator -= (const Vector& A) {
        coords[0] -= A[0];
        coords[1] -= A[1];
        coords[2] -= A[2];
        return *this;
    }
    Vector& operator -= (const double& a) {
        coords[0] -= a;
        coords[1] -= a;
        coords[2] -= a;
        return *this;
    }
    Vector& operator *= (const double& a) {
        coords[0] *= a;
        coords[1] *= a;
        coords[2] *= a;
        return *this;
    }
    Vector& operator /= (const double& a) {
        coords[0] /= a;
        coords[1] /= a;
        coords[2] /= a;
        return *this;
    }

    double norm2() {
        return sqr(coords[0]) + sqr(coords[1]) + sqr(coords[2]);
    }
    Vector normalize() {
        double norm = sqrt(this->norm2());
        coords[0] /= norm;
        coords[1] /= norm;
        coords[2] /= norm;
        return *this;
    }

private:
    double coords[3];
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
Vector operator*(const Vector &A, double b) {
    return Vector(A[0] * b, A[1] * b, A[2] * b);
}
Vector operator*(double b, const Vector &A) {
    return Vector(A[0] * b, A[1] * b, A[2] * b);
}
Vector operator/(const Vector &A, double b) {
    return Vector(A[0] / b, A[1] / b, A[2] / b);
}
Vector operator/(double &b, const Vector &A) {
    return Vector(b / A[0], b / A[1], b / A[2]);
}
Vector operator/(const Vector &A, const Vector& B) {
    return Vector(A[0] / B[0], A[1] / B[1], A[2] / B[2]);
}
double dot(const Vector& A, const Vector& B) {
    return A[0] * B[0] + A[1] * B[1] + A[2] * B[2];
};
Vector cross(const Vector& A, const Vector& B) {
    return Vector(A[1] * B[2] - B[1] * A[2], A[2] * B[0] - B[2] * A[0], A[0] * B[1] - B[0] * A[1]);
}

Vector randomCos(const Vector& N) {
    int tid = omp_get_thread_num();
    double r1 = uniform(engine[tid]);
    double r2 = uniform(engine[tid]);
    double r = sqrt(1-r2);
    double x = cos(2.*M_PI*r1) * r;
    double y = sin(2.*M_PI*r1) * r;
    double z = sqrt(r2);

    Vector T1;
    // Removing smallest compound to avoid round errors (+ limited impact)
    if (abs(N[0]) <= abs(N[1]) && abs(N[0]) <= abs(N[2])) {
        T1 = Vector(0, -N[2], N[1]);
    }
    else {
        if (abs(N[1]) <= abs(N[2]) && abs(N[1]) <= abs(N[0])) {
            T1 = Vector(-N[2], 0, N[0]);
        }
        else {
            T1 = Vector(-N[1], N[0], 0);
        }
    }
    T1.normalize();
    Vector T2 = cross(N, T1);

    return x * T1 + y * T2 + z * N;
}

Vector max(const Vector& A, const Vector& B) {
    return Vector(std::max(A[0], B[0]), std::max(A[1], B[1]), std::max(A[2], B[2]));
}

double max(const Vector& A) {
    return std::max(A[0], std::max(A[1], A[2]));
}

Vector min(const Vector& A, const Vector& B) {
    return Vector(std::min(A[0], B[0]), std::min(A[1], B[1]), std::min(A[2], B[2]));
}

double min(const Vector& A) {
    return std::min(A[0], std::min(A[1], A[2]));
}

int idMax(const Vector& A) {
    if (A[0] >= A[1]) {
        if (A[0] >= A[2]) {
            return 0;
        }
        else {
            return 2;
        }
    }
    else {
        if (A[1] >= A[2]) {
            return 1;
        }
        else {
            return 2;
        }
    }
}

int idMin(const Vector& A) {
    if (A[0] <= A[1]) {
        if (A[0] <= A[2]) {
            return 0;
        }
        else {
            return 2;
        }
    }
    else {
        if (A[1] <= A[2]) {
            return 1;
        }
        else {
            return 2;
        }
    }
}

#endif