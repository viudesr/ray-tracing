#ifndef OBJECT_HPP
#define OBJECT_HPP
#include "ray.hpp"

class Object {
public:
    Object(const Vector& rho, bool mirror, bool transparent, double n, bool light, double lightIntensity) : rho(rho), mirror(mirror), transparent(transparent), n(n), light(light), lightIntensity(lightIntensity) {};

    virtual bool intersect(const Ray& ray, double& t, Vector& N, Vector& P) const=0;
    Vector rho;
    bool mirror, transparent, light;
    double n, lightIntensity;
};

class Sphere : public Object {
public:
    Sphere(const Vector& O, double r, const Vector& rho, bool mirror = false, bool transparent = false, double n = 1.4, bool light = false, double lightIntensity = 0.) : Object(rho, mirror, transparent, n, light, lightIntensity), radius(r), origin(O) {}

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
	double radius;
	Vector origin;
};

class Triangle: public Object {
public:
	Triangle(const Vector& A, const Vector& B, const Vector& C, const Vector& rho, bool mirror=false, bool transparent=false, double n=1.4, bool light=false, double lightIntensity = 0.) : Object(rho, mirror, transparent, n, light, lightIntensity), A(A), B(B), C(C) {};

	bool intersect(const Ray& r, double &t, Vector &P, Vector &N) const {
		N = cross(B - A, C - A);
		N.normalize();
		t = dot(C - r.origin, N) / dot(r.dir, N);
		if (t < 0) return false;

		P = r.origin + t * r.dir;
		Vector u = B - A;
		Vector v = C - A;
		Vector w = P - A;
		double m11 = u.norm2();
		double m12 = dot(v,u);
		double m22 =  v.norm2();
		double detm = m11*m22 - sqr(m12);

		double b11 = dot(u,w);
		double b21 = dot(w,v);
		double detb = b11*m22 - b21*m12;
		double beta = detb / detm;

		double g12 = b11;
		double g22 = b21;
		double detg = m11 * g22 - m12 * g12;
		double gamma = detg / detm;

		double alpha = 1 - beta - gamma;
		if (alpha < 0 || alpha > 1) return false;
		if (beta < 0 || beta > 1) return false;
		if (gamma < 0 || gamma > 1) return false;

		return true;
	}
	Vector A, B, C;
};

#endif