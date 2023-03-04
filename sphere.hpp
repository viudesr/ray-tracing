#ifndef SPHERE_HPP
#define SPHERE_HPP
#include "vector.hpp"

class Sphere {
public:
    Sphere(Vector& O, double &r) {
        origin = O;
        radius = r;
    }
    Vector origin;
    double radius;
};

#endif