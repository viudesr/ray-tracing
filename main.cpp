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
    int W = 1024;
    int H = 1024;
    std::vector<unsigned char> image(W * H * 3, 0);

    double fov = 55 * M_PI / 180.;
    double gamma = 2.2;
    Vector camera(0., 0., 55.);
    Scene scene;
    Sphere S(Vector(0.,0.,0.), 10.);
    Sphere S2(Vector(20.,0.,0.), 10.);

    scene.addSphere(S);
    scene.addSphere(S2);

#pragma omp parallel for
    for (int i = 0; i < H; i++) {
        for (int j = 0; j < W; j++) {

            int x = j - W / 2; // x coordinate of pixel
            int y = i - H / 2; // y coordinate of pixel
            int z = -W / (2 * tan(fov / 2)); // z coordinate of pixel

            Vector u(x, y, z);
            u.normalize();
            Vector color(0.,0.,0.);
            if (scene.intersect(Ray(camera, u))) {
                color = Vector(255.,255.,255.);
            }

            image[(i*W + j) * 3 + 0] = color[0];   // RED
            image[(i*W + j) * 3 + 1] = color[1];  // GREEN
            image[(i*W + j) * 3 + 2] = color[2];  // BLUE

        }
    }
    stbi_write_png("figures/exemple.png", W, H, 3, &image[0], 0);

    return 0;
}