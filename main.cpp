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
    int Nrays = 32;
    std::vector<unsigned char> image(W * H * 3, 0);

    double fov = 55 * M_PI / 180.;
    double gamma = 2.2;
    Vector camera(0., 0., 55.);
    Scene scene;
    Vector lightSource1(20,25,20);
    Vector lightSource2(-20, 25, -50);
    Sphere Slum(lightSource1, 5., Vector(1.,1.,1.), false, false, 1.4, true, 3E8);
    Sphere Slum2(lightSource2, 5., Vector(1.,1.,1.), false, false, 1.4, true, 3E8);
    Sphere S(Vector(0.,0.,-55.), 10., Vector(0.9,0.5,0.5));
    Sphere S2(Vector(20.,0.,-55.), 10., Vector(0.9,0.5,0.5), true);
    Sphere S3(Vector(-20.,0.,-55.), 10., Vector(0.9,0.5,0.5), false, true);
    Sphere S_bottom(Vector(0.,-1000.,0.), 990., Vector(0.6,0.3,0.6));
    Sphere S_top(Vector(0.,1000.,0.), 970., Vector(0.1,0.3,0.9));
    Sphere S_left(Vector(-1000.,0.,0.), 970., Vector(0.5,0.5,0.6));
    Sphere S_right(Vector(1000.,0.,0.), 970., Vector(0.9,0.3,0.1));
    Sphere S_back(Vector(0.,0.,-1000.), 920., Vector(0.7,0.4,0.6));
    Sphere S_front(Vector(0.,0.,1000), 940., Vector(0.1,0.1,0.8));

    scene.addSphere({Slum, Slum2});
    scene.addSphere({S, S2, S3, S_bottom, S_top, S_left, S_right, S_back, S_front});

#pragma omp parallel for
    for (int i = 0; i < H; i++) {
        int tid = omp_get_thread_num();
        for (int j = 0; j < W; j++) {

            Vector color(0.,0.,0.);
            for (int k = 0; k < Nrays; k++) {
                // Anti-aliasing (averaged small random variation of direction)
                double r1 = uniform(engine[tid]);
                double r2 = uniform(engine[tid]);
                double r = sqrt(-2*log(r1));
                double gx = r * cos(2 * M_PI * r2) * 0.7;
                double gy = r * sin(2 * M_PI * r2) * 0.7;

                // Ray coordinates
                double x = j - W / 2 + gx; // x coordinate of pixel
                double y = i - H / 2 + gy; // y coordinate of pixel
                double z = -W / (2. * tan(fov / 2)); // z coordinate of pixel

                Vector u(x, y, z);
                u.normalize();

                // Shooting ray(s)
                color += scene.getColor(Ray(camera,u), 5, false);
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