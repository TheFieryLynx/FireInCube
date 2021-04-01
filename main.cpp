#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>
#include "geom.h"
#include "functions.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

Pixel cast_ray(const Ray& ray, const Cube& cube)
{
    float t;
    if (!cube.ray_intersect(ray, t)) {
        //return Pixel(255 * (ray.dir.x / 2 + 0.5), 255 * (ray.dir.y / 2 + 0.5), 255 * (ray.dir.z / 2 + 0.5));
        return Pixel(0, 200, 255);
    }
    return Pixel (140, 50, 50);
}

void render(const Cube& cube, const Vec3D& camera) {
    const size_t width = 768;
    const size_t height = 768;
    Pixel* imageBuffer = new Pixel[width * height];
    for (size_t j = 0; j < height; ++j) {
        for (size_t i = 0; i < width; ++i) {
            float x = -(float)width / 2 + i + 0.5;
            float y = -(float)height / 2 + j + 0.5;
            float z = 39;
            Vec3D current_point(x, y, z);
            Vec3D dir = Vec3D(current_point.x - camera.x, current_point.y - camera.y, current_point.z - camera.z).normalize();
            Ray ray = Ray(camera, dir);
            imageBuffer[i + j * width] = cast_ray(ray, cube);
        
        }
    }
    stbi_flip_vertically_on_write(0);
    stbi_write_png("../out.png", width, height, sizeof(Pixel), imageBuffer, sizeof(Pixel) * width);
}

int main()
{
    Cube cube(Vec3D(-100, -100, 40), Vec3D(100, 100, 60));
    Vec3D camera(-300, 300, 0);
    render(cube, camera);
    return 0;
}