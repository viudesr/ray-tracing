#ifndef OBJECT_HPP
#define OBJECT_HPP
#include "ray.hpp"

class Sphere {
public:
    Sphere(Vector& O, double& r) {
        origin = O;
        radius = r;
    }

    bool intersect(const Ray& ray, double& t) {
        double b = dot(ray.dir, ray.origin - this->origin);
        double delta =  sqr(b) - ((ray.origin - this->origin).norm2() - sqr(this->radius));
        if (delta < 0) {
            return false;
        }
        double t0 = 2 * (- b - sqrt(delta));
        double t1 = 2 * (- b + sqrt(delta));

        if (t1 < 0) {
            return false;
        }

        if (t0 > 0) {
            t = t0;
        }
        else {
            t = t1;
        }

        return true;
    }
    Vector origin;
    double radius;
};

#endif