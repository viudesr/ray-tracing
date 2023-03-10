#ifndef RAY_HPP
#define RAY_HPP
#include "vector.hpp"

class Ray {
public:
    Ray(const Vector& Origin, const Vector& u) : origin(Origin), dir(u) {}
    Vector origin;
    Vector dir;
};
#endif