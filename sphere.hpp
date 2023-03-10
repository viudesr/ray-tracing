#ifndef SPHERE_HPP
#define SPHERE_HPP
#include "vector.hpp"

class Sphere {
public:
    Sphere(Vector& O, double& r) {
        origin = O;
        radius = r;
    }

    bool intersect(const Ray& ray, double& t) {
        double delta =  sqr(2 * dot(ray.dir, ray.origin - this->origin)) - 4 * ((ray.origin - this->origin).norm2() - sqr(this->radius));
        if (delta > 0) {
            t = - dot(ray.dir, ray.origin - this->origin) - sqrt(delta) / 2;
            return true;
        }
    }
    Vector origin;
    double radius;
};

#endif