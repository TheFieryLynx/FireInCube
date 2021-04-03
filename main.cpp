#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>
#include "geom.h"


#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

float cosinus(const Vec3D& a, const Vec3D& b) 
{
    float m = dot(a, b);
    float n1 = std::sqrt(std::pow(a.x, 2) + std::pow(a.y, 2) + std::pow(a.z, 2));
    float n2 = std::sqrt(std::pow(b.x, 2) + std::pow(b.y, 2) + std::pow(b.z, 2));
    return m / (n1 * n2);
}

bool scene_sphere_intersect(const Ray& ray, const Sphere &sphere, float &dist_i, Vec3D &N, Material &material) {
    float sphere_dist = std::numeric_limits<float>::max();
    if (sphere.ray_intersect(ray.orig, ray.dir, dist_i) && dist_i < sphere_dist) {
        sphere_dist = dist_i;
        Vec3D hit = ray.orig + ray.dir * dist_i;
        N = (hit - sphere.center).normalize();
        material = sphere.material;
        return true;
    }
    return false;
}

float scene_cube_intersect(const Ray& ray, const Cube &cube, float& t, Vec3D &N, Material &material, const std::vector<Light>& lights) {
    float cube_dist = std::numeric_limits<float>::max();
    if (cube.ray_intersect(ray, t)) {
        cube_dist = t;
        Vec3D point = ray.orig + ray.dir * t;
        //std::cout << point << std::endl;
        //N = (point - cube.center).normalize();
        N = cube.normInPoint(point);
        //std::cout << N << std::endl;
        float diffuse_light_intensity = 0.2;
        for (size_t i = 0; i < lights.size(); ++i) {
            Vec3D light_dir = (lights[i].position - point);
            //std::cout << point << std::endl;
            diffuse_light_intensity  += lights[i].intensity * std::max(0.f, cosinus(N, light_dir));  //fabs(cosinus(light_dir, N)); //std::max(0.f, light_dir * N);
        }
        material = cube.material;
        material.diffuse_color = material.diffuse_color * diffuse_light_intensity;
        // Pixel p = Pixel(N.x + 1, N.y + 1, N.z + 1) * 255;
        // material.diffuse_color = p;
        return true;
    }
    return false;
    
    //std::cout << diffuse_light_intensity << std::endl;
}

Pixel cast_ray(const Ray& ray, const std::vector<Sphere>& spheres, const Cube& cube, const std::vector<Light>& lights)
{
    float t_sphere = std::numeric_limits<float>::max();
    Vec3D sphere_N;
    Material sphere_material;
    int sphere_flag = 1, cube_flag = 0;
    if (!scene_sphere_intersect(ray, spheres[0], t_sphere, sphere_N, sphere_material)) {
        sphere_flag = 1;
        //return Pixel(0, 200, 255); // background color
    }
    float t_cube = std::numeric_limits<float>::max();
    Vec3D cube_N;
    Material cube_material;
    if (!scene_cube_intersect(ray, cube, t_cube, cube_N, cube_material, lights)) {
        cube_flag = 1;
        //return Pixel(0, 200, 255); // background color
    }
    // float d_sphere = std::numeric_limits<float>::max(), d_cube = std::numeric_limits<float>::max();
    // d_sphere = std::sqrt(std::pow(sphere_point.x - ray.orig.x, 2) + std::pow(sphere_point.y - ray.orig.y, 2) + std::pow(sphere_point.z - ray.orig.z, 2));
    // d_cube = std::sqrt(std::pow(cube_point.x - ray.orig.x, 2) + std::pow(cube_point.y - ray.orig.y, 2) + std::pow(cube_point.z - ray.orig.z, 2));

    //std::cout << d_sphere << " " << d_cube << std::endl;

    if (cube_flag && sphere_flag) {
        return Pixel(0, 200, 255);
    }
    
    //std::cout << t_cube << " " << t_sphere << std::endl;

    //if (t_cube < t_sphere) {
        //return Pixel(78, 26, 26);
        return cube_material.diffuse_color;
    //}
    //return Pixel(102, 102, 78);
    return sphere_material.diffuse_color;

    //return Pixel(140, 50, 50);
    
}

void render(const std::vector<Sphere>& spheres, const Cube& cube, const Vec3D& camera, const std::vector<Light>& lights) {
    const size_t width = 768;
    const size_t height = 768;
    Pixel* imageBuffer = new Pixel[width * height];
    for (size_t j = 0; j < height; ++j) {
        for (size_t i = 0; i < width; ++i) {
            // float x = (float)i / width + 0.5f;
            // float y = (float)j / height + 0.5f;
            // float z = 3;
            float x = -(float)width / 2 + i + 0.5;
            float y = -(float)height / 2 + j + 0.5;
            float z = 768;
            Vec3D current_point(x, y, z);
            Vec3D dir = Vec3D((current_point.x - camera.x), (current_point.y - camera.y), (current_point.z - camera.z)).normalize();
            Ray ray = Ray(camera, dir);
            imageBuffer[i + j * width] = cast_ray(ray, spheres, cube, lights);
        
        }
    }
    stbi_flip_vertically_on_write(0);
    stbi_write_png("../out.png", width, height, sizeof(Pixel), imageBuffer, sizeof(Pixel) * width);
}

int main()
{
    Material ivory(Pixel(102, 102, 78));
    Material red_rubber(Pixel(178, 26, 26));

    std::vector<Sphere> spheres;

    Cube cube(Vec3D(0, 0, 750), 100, red_rubber); 

    spheres.push_back(Sphere(Vec3D(0, -200, 650), 50, ivory));

    //Cube cube(Vec3D(-100, -100, 50), Vec3D(100, 100, 250));
    Vec3D camera(-300, -500, -300);

    std::vector<Light> lights;
    lights.push_back(Light(Vec3D(-200, -250, 850), 1));

    render(spheres, cube, camera, lights);
    return 0;
}