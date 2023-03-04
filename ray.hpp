#ifndef RAY_HPP
#define RAY_HPP
#include "vector.hpp"

class Ray {
public:
    Ray(double& Origin, Vector& u) {
        origin = Origin;
        u = u;
    }
    double origin;
    Vector u;
};
#endif