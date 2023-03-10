#ifndef SCENE_HPP
#define SCENE_HPP

class Scene {
public:
    Scene() {};
    void addSphere(const Sphere& s) {
        objects.push_back(s);
    };

    bool intersect(const Ray& ray) {
        double t;
        if (objects[0].intersect(ray, t)) {
            return true;
        }
        else {
            return false;
        }
        for (int i=0; i < objects.size(); i++) {
            double t;
            if (objects[i].intersect(ray, t)) {
                return true;
            }
        }
        return false;
    }

    std::vector<Sphere> objects;
};

#endif