#ifndef OBJECT_HPP
#define OBJECT_HPP
#include "ray.hpp"

class Object {
public:
    Object() {};

    virtual bool intersect(const Ray& ray, double& t) const=0;
};

class Sphere : public Object {
public:
    Sphere(const Vector& O, double r) : origin(O), radius(r) {}

    bool intersect(const Ray& ray, double& t) const {
        double b = 2 * dot(ray.dir, ray.origin - this->origin);
        double delta =  sqr(b) - 4 * ((ray.origin - this->origin).norm2() - sqr(this->radius));
        if (delta < 0) {
            return false;
        }
        double t0 = (- b - sqrt(delta)) / 2;
        double t1 = (- b + sqrt(delta)) / 2;

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