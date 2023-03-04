#ifndef RAY_HPP
#define RAY_HPP
#include "vector.hpp"

class Ray {
public:
    Ray(Vector& Origin, Vector& u) {
        origin = Origin;
        u = u;
    }
    Vector origin;
    Vector u;
};
#endif