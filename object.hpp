#ifndef OBJECT_HPP
#define OBJECT_HPP
#include "ray.hpp"

class Object {
public:
    Object(const Vector& rho, bool mirror, bool transparent, double n) : rho(rho), mirror(mirror), transparent(transparent), n(n) {};

    virtual bool intersect(const Ray& ray, double& t, Vector& N, Vector& P) const=0;
    Vector rho;
    bool mirror, transparent;
    double n;
};

class Sphere : public Object {
public:
    Sphere(const Vector& O, double r, const Vector& rho, bool mirror = false, bool transparent = false, double n = 1.4) : Object(rho, mirror, transparent, n), origin(O), radius(r) {}

    bool intersect(const Ray& ray, double& t, Vector& N, Vector& P) const {
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

        P = ray.origin + t * ray.dir;
        N = P - this->origin;
        N.normalize();

        return true;
    }
    Vector origin;
    double radius;
};

#endif