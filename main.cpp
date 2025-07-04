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
    int W = 512;
    int H = 512;
    int Nrays = 32;
    std::vector<unsigned char> image(W * H * 3, 0);

    // Camera settings
    double fov = 55 * M_PI / 180.;
    double gamma = 2.2;
    double aperture = 1.;
    double focal_length = 95.;
    Vector camera(0., 0., 55.);

    // Scene setup
    Scene scene;
    Vector lightSource1(20,25,20);
    Vector lightSource2(-20, 25, -50);
    Sphere Slum(lightSource1, 5., Vector(1.,1.,1.), false, false, 1.4, true, 3E8);
    Sphere Slum2(lightSource2, 5., Vector(1.,1.,1.), false, false, 1.4, true, 3E8);
    Sphere S(Vector(0.,0.,-55.), 10., Vector(0.9,0.5,0.5));
    Sphere S2(Vector(20.,0.,-55.), 10., Vector(0.9,0.5,0.5), true);
    Sphere S3(Vector(-20.,0.,-55.), 10., Vector(0.9,0.5,0.5), false, true);
    Sphere S4(Vector(0., 0., -20), 10., Vector(1.,1.,1.), false, true);
    Sphere S_bottom(Vector(0.,-1000.,0.), 990., Vector(0.6,0.3,0.6));
    Sphere S_top(Vector(0.,1000.,0.), 970., Vector(0.1,0.3,0.9));
    Sphere S_left(Vector(-1000.,0.,0.), 970., Vector(0.5,0.5,0.6));
    Sphere S_right(Vector(1000.,0.,0.), 970., Vector(0.9,0.3,0.1));
    Sphere S_back(Vector(0.,0.,-1000.), 920., Vector(0.7,0.4,0.6));
    Sphere S_front(Vector(0.,0.,1000), 940., Vector(0.1,0.1,0.8));

    TriangleMesh mesh(Vector(1.,0.,0.));
    mesh.readOBJ("models/cat.obj");
    mesh.transform(0.5, Vector(0.,-10.,10.));
    mesh.computeBVH();

    Triangle T1(Vector(-10,-10,-20), Vector(50,-10,-20), Vector(0,50,-19), Vector(1.,0.,0.));
    //mesh.makeTri(T1);
    //mesh.computeBVH();
    
    // Adding objects to scene
    scene.addSphere({Slum, Slum2, S, S2, S3, S_bottom, S_top, S_left, S_right, S_back, S_front});
    scene.addMesh(mesh);
    //scene.addTriangle(T1);

#pragma omp parallel for
    for (int i = 0; i < H; i++) {
        int tid = omp_get_thread_num();
        for (int j = 0; j < W; j++) {

            // Ray coordinates
            double x = j - W / 2 + 0.5; // x coordinate of pixel
            double y = i - H / 2 + 0.5; // y coordinate of pixel
            double z = -W / (2. * tan(fov / 2)); // z coordinate of pixel
            Vector color(0.,0.,0.);

            // Basic ray
            Ray ray(camera, Vector(x, y, z));

            for (int k = 0; k < Nrays; k++) {
                // Anti-aliasing (averaged small random variation of direction)
                antiAliasing(ray, camera, x, y, z);

                // Camera model
                //cameraModel(ray, focal_length, aperture);

                // Shooting ray(s)
                color += scene.getColor(ray, 5, false);
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