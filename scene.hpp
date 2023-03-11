#ifndef SCENE_HPP
#define SCENE_HPP

#include <initializer_list>



class Scene {
public:
    Scene() {
        I = 1E6;
        lightSource = Vector(-10,20,40);
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

    Vector getColor(const Ray& ray) {
        double t;
        Vector N, P;
        int id;

        Vector color(0., 0., 0.);
        if (this->intersect(ray, t, N, P, id)) {
            Vector object2Light = lightSource - P;
            double dist2Light = object2Light.norm2();
            object2Light.normalize();
            color = I * objects[id].rho / M_PI * std::max(0., dot(object2Light, N + 0.001)) / dist2Light;
        }
        return color;
    }

    double I;
    Vector lightSource;
    std::vector<Sphere> objects;
};

#endif