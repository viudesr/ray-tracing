#ifndef RAY_HPP
#define RAY_HPP
#include "vector.hpp"

class Ray {
public:
    Ray(const Vector& Origin, const Vector& u) : origin(Origin), dir(u) {}
    Vector origin;
    Vector dir;
};

static inline Ray reflectedRay(const Ray& ray, const Vector& newOrigin, const Vector& N) {
    Vector rayN = dot(ray.dir, N) * N;
    return Ray(newOrigin + 0.001 * N, ray.dir - 2 * rayN);
}

static inline Ray refractedRay(const Ray& ray, const Vector& newOrigin, const Vector& N, double n1, double n2) {
    double rayN = dot(ray.dir, N);
    Vector Tt = n1 / n2 * (ray.dir - rayN * N);
    Vector Tn = sqrt(1 - sqr(n1 / n2) * (1 - sqr(rayN))) * - N;
    return Ray(newOrigin - 0.001 * N, Tt + Tn);
}

static inline Ray randomRay(const Vector& newOrigin, Vector& N) {
    return Ray(newOrigin + 0.001 * N, randomCos(N));
}

static inline Vector boxMuller(double std) {
    int tid = omp_get_thread_num();
    double r1 = uniform(engine[tid]);
    double r2 = uniform(engine[tid]);
    double r = sqrt(-2 * log(r1));
    double x = r * cos(2 * M_PI * r2) * std;
    double y = r * sin(2 * M_PI * r2) * std;
    double z = sqrt(r2);
    return Vector(x, y, z);
}

static inline Ray cameraModel(Ray& ray, double focal_length, double aperture) {
    Vector dx = boxMuller(aperture);
    Vector focusPoint = ray.origin + focal_length * ray.dir;
    Vector newOrigin = ray.origin + Vector(dx[0], dx[1], 0.);
    Vector newDir = focusPoint - newOrigin;
    newDir.normalize();
    return Ray(newOrigin, newDir);
}
#endif