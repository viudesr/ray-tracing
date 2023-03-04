
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include "vector.hpp"
#include "ray.hpp"
#include "sphere.hpp"
#include "scene.hpp"

int main() {
    int W = 1024;
    int H = 1024;
    double fov = 55;
    for (int i=0; i<H; i++) {
        for (int j=0; j<W; j++) {
            int x = j - W/2;
            int y = i - H/2;
            int z = -W/(2*tan(fov/2));
        }
    }
    return 0;
}