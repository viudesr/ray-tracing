#define _CRT_SECURE_NO_WARNINGS 1
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include <vector>
#include <omp.h>
#include <iostream>

#include "vector.hpp"
#include "ray.hpp"
#include "object.hpp"
#include "scene.hpp"

int main() {
    int W = 2048;
    int H = 2048;
    int Nrays = 32;
    std::vector<unsigned char> image(W * H * 3, 0);

    double fov = 55 * M_PI / 180.;
    double gamma = 2.2;
    Vector camera(0., 0., 55.);
    Scene scene;
    Sphere S(Vector(0.,0.,-55.), 10., Vector(0.9,0.5,0.5));
    Sphere S2(Vector(20.,0.,-55.), 10., Vector(0.9,0.5,0.5), true);
    Sphere S3(Vector(-20.,0.,-55.), 10., Vector(0.9,0.5,0.5), false, true);
    Sphere S_bottom(Vector(0.,-1000.,0.), 990., Vector(0.6,0.3,0.6));
    Sphere S_top(Vector(0.,1000.,0.), 970., Vector(0.1,0.3,0.9));
    Sphere S_left(Vector(-1000.,0.,0.), 970., Vector(0.5,0.5,0.6));
    Sphere S_right(Vector(1000.,0.,0.), 970., Vector(0.9,0.3,0.1));
    Sphere S_back(Vector(0.,0.,-1000.), 920., Vector(0.7,0.4,0.6));
    Sphere S_front(Vector(0.,0.,1000), 940., Vector(0.1,0.1,0.8));

    scene.addSphere({S, S2, S3, S_bottom, S_top, S_left, S_right, S_back, S_front});

#pragma omp parallel for
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {

            int x = j - W / 2; // x coordinate of pixel
            int y = i - H / 2; // y coordinate of pixel
            int z = -W / (2 * tan(fov / 2)); // z coordinate of pixel

            Vector u(x, y, z);
            u.normalize();

            Vector color(0.,0.,0.);
            for (int k = 0; k < Nrays; k++) {
                color += scene.getColor(Ray(camera,u), 5);
            }
            color /= Nrays;

            image[((H - i - 1) * W + j) * 3 + 0] = std::min(pow(color[0], 1./gamma), 255.);   // RED
            image[((H - i - 1) * W + j) * 3 + 1] = std::min(pow(color[1], 1./gamma), 255.);  // GREEN
            image[((H - i - 1) * W + j) * 3 + 2] = std::min(pow(color[2], 1./gamma), 255.);  // BLUE

        }
    }
    stbi_write_png("figures/exemple.png", W, H, 3, &image[0], 0);

    return 0;
}