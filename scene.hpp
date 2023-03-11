#ifndef SCENE_HPP
#define SCENE_HPP

#include <initializer_list>

class Scene {
public:
    Scene() {
        I = 1E9;
        lightSource = Vector(20,30,50);
    };

    template <class Sphere>
    void addSphere(std::initializer_list <Sphere> list) {
        for (const Sphere& s : list) {
            objects.push_back(s);
        }
    };

    void addSphere(const Sphere& s) {
        objects.push_back(s);
    }

    bool intersect(const Ray& ray, double& t, Vector& N, Vector& P, int& id) {
        /* Checks which object has the closest intersection with this ray, and returns intersection parameters*/
        bool intersection = false;
        
        double local_t;
        t = std::numeric_limits<double>::max();

        Vector local_N, local_P;
        for (int i=0; i < objects.size(); i++) {
            if (objects[i].intersect(ray, local_t, local_N, local_P)) {
                intersection = true;
                if (local_t < t) {
                    t = local_t;
                    N = local_N;
                    P = local_P;
                    id = i;
                }
            }
        }
        return intersection;
    }

    Vector getColor(const Ray& ray, int n_bounces) {
        /* Returns color vector associated to light for this ray */

        Vector color(0., 0., 0.);

        // Case of end of bounces
        if (n_bounces <= 0) {
            return color;
        }
        double t;
        Vector N, P;
        int id;

        if (this->intersect(ray, t, N, P, id)) {
            // Case of intersection between this ray and an object

            //Handling mirror case
            if (objects[id].mirror) {
                return getColor(reflectedRay(ray, P, N), n_bounces - 1);
            }

            //Handling transparency
            if (objects[id].transparent) {
                double n1 = 1.;
                double n2 = objects[id].n;
                Vector N_temp = N;
                if (dot(ray.dir, N) > 0) {
                    // case where ray is leaving the sphere
                    std::swap(n1, n2);
                    N_temp *= -1;
                }
                double cos2_2 = 1 - n1 / n2 * (1 - sqr(dot(ray.dir, N_temp)));
                if (cos2_2 >= 0) {
                    return getColor(refractedRay(ray, P, N_temp, n1, n2), n_bounces - 1);
                }
                else {
                    return getColor(reflectedRay(ray, P, N), n_bounces - 1);
                }
            }

            // Direct lighting calculation
            
            Vector object2Light = lightSource - P;
            double dist2_2Light = object2Light.norm2();
            object2Light.normalize();

            // Visibility term computation

            bool visibility = true;
            double tprime;
            Vector Pprime, Nprime;
            int idprime;

            bool intersectSecondaryObject = intersect(Ray(P + 0.001 * N, object2Light), tprime, Pprime, Nprime, idprime);
            if (intersectSecondaryObject && (sqr(tprime) < dist2_2Light)) {
                visibility = false;
            }

            // Indirect lighting calculation

            Vector indirect = objects[id].rho * getColor(randomRay(P, N), n_bounces - 1);

            // color affectation

            color = visibility * I * objects[id].rho / M_PI * std::max(0., dot(object2Light, N + 0.001)) / dist2_2Light;
            color += indirect;
        }
        return color;
    }

    double I;
    Vector lightSource;
    std::vector<Sphere> objects;
};

#endif