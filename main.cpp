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

float sinus(const Vec3D& a, const Vec3D& b) 
{
    float cos = cosinus(a, b);
    return std::sqrt(1 - std::pow(cos, 2));
}

bool scene_sphere_intersect(const Ray& ray, const Sphere &sphere, float &dist_i, Vec3D &N, Material &material, const std::vector<Light>& lights) {
    if (sphere.ray_intersect(ray.orig, ray.dir, dist_i)) {
        Vec3D point = ray.orig + ray.dir * dist_i;
        N = (point - sphere.center).normalize();
        float diffuse_light_intensity = 0.2;
        
        for (int i = 0; i < int(lights.size()); ++i) {
            Vec3D light_dir = (lights[i].position - point).normalize();
            diffuse_light_intensity  += lights[i].intensity * std::max(0.f, cosinus(N, light_dir));  //fabs(cosinus(light_dir, N)); //std::max(0.f, light_dir * N);
        }
        material = sphere.material;
        material.diffuse_color = material.diffuse_color * diffuse_light_intensity;
        return true;
    } 
    return false;
}

float scene_cube_intersect(const Ray& ray, const Cube &cube, float& t, Vec3D &N, Material &material, const std::vector<Light>& lights) {
    if (cube.ray_intersect(ray, t)) {
        Vec3D point = ray.orig + ray.dir * t;
        N = cube.normInPoint(point);
        float diffuse_light_intensity = 0.2;
        for (int i = 0; i < int(lights.size()); ++i) {
            Vec3D light_dir = (lights[i].position - point).normalize();
            diffuse_light_intensity  += lights[i].intensity * std::max(0.f, cosinus(N, light_dir));  //fabs(cosinus(light_dir, N)); //std::max(0.f, light_dir * N);
        }
        material = cube.material;
        material.diffuse_color = material.diffuse_color * diffuse_light_intensity;
        return true;
    }
    return false;
}

Pixel cast_ray(const Ray& ray, const std::vector<Sphere>& spheres, const Cube& cube, const std::vector<Light>& lights)
{
    std::vector<int> flags = { 1, 1, 1 };
    std::vector<float> t_s;
    std::vector<Material> materials;

    float t_cube = std::numeric_limits<float>::max();
    Vec3D cube_N;
    Material cube_material;
    int cube_flag = 1;
    if (!scene_cube_intersect(ray, cube, t_cube, cube_N, cube_material, lights)) {
        flags[0] = 0;
        t_cube = -1.0;
        cube_material = Material(Pixel(0, 200, 255));
        //return Pixel(0, 200, 255); // background color
    }
    materials.push_back(cube_material);
    t_s.push_back(t_cube);

    float t_sphere1 = std::numeric_limits<float>::max();
    Vec3D sphere1_N;
    Material sphere1_material;
    if (!scene_sphere_intersect(ray, spheres[0], t_sphere1, sphere1_N, sphere1_material, lights)) {
        flags[1] = 0;
        t_sphere1 = -1.0;
        sphere1_material = Material(Pixel(0, 200, 255));
        //return Pixel(0, 200, 255); // background color
    }
    materials.push_back(sphere1_material);
    t_s.push_back(t_sphere1);

    float t_sphere2 = std::numeric_limits<float>::max();
    Vec3D sphere2_N;
    Material sphere2_material;
    if (!scene_sphere_intersect(ray, spheres[1], t_sphere2, sphere2_N, sphere2_material, lights)) {
        flags[2] = 0;
        t_sphere2 = -1.0;
        sphere2_material = Material(Pixel(0, 200, 255));
    }
    materials.push_back(sphere2_material);
    t_s.push_back(t_sphere2);
    
    float tmp = std::numeric_limits<float>::max();
    int i_tmp = -1;
    for (int i = 0; i < flags.size(); ++i) {
        if (flags[i] != 0) {
            if (t_s[i] < tmp) {
                tmp = t_s[i];
                i_tmp = i;
            }
        }
    }
    if (i_tmp == -1) {
        return Pixel(0, 200, 255);
    }
    return materials[i_tmp].diffuse_color;

    // if (!cube_flag && !sphere_flags[0] && !sphere_flags[1]) {
    //     return Pixel(0, 200, 255);
    // }
    
    //std::cout << t_cube << " " << t_sphere << std::endl;


    // if (cube_flag && (!sphere_flag || (t_cube < t_sphere))) {
    //     //return Pixel(78, 26, 26);
    //     return cube_material.diffuse_color;
    // }
    // return sphere_material.diffuse_color;

    //return Pixel(140, 50, 50);
    
}

Vec3D Rotation_x(Vec3D point, float a) {
    float EPS = 0.00001;
    float k = M_PI / 180;
    float cos_x = std::cos(k * a);
    float sin_x = std::sin(k * a);
    Vec3D new_point;
    new_point.x = point.x;
    if (fabs(new_point.x) < EPS) {
        new_point.x = 0.;
    }
    new_point.y = point.y * cos_x + point.z * sin_x;
    if (fabs(new_point.y) < EPS) {
        new_point.y = 0.;
    }
    new_point.z = - point.y * sin_x + point.z * cos_x;
    if (fabs(new_point.z) < EPS) {
        new_point.z = 0.;
    }
    return new_point;
}

Vec3D Rotation_y(Vec3D point, float a) {
    float EPS = 0.00001;
    float k = M_PI / 180;
    float cos_y = std::cos(k * a);
    float sin_y = std::sin(k * a);
    Vec3D new_point;
    new_point.x = point.x * cos_y + point.z * sin_y;
    if (fabs(new_point.x) < EPS) {
        new_point.x = 0.;
    }
    new_point.y = point.y;
    if (fabs(new_point.y) < EPS) {
        new_point.y = 0.;
    }
    new_point.z = -point.x * sin_y + point.z * cos_y;
    if (fabs(new_point.z) < EPS) {
        new_point.z = 0.;
    }
    return new_point;
}

Vec3D Rotation_z(Vec3D point, float a) {
    float EPS = 0.00001;
    float k = M_PI / 180;
    float cos_z = std::cos(k * a);
    float sin_z = std::sin(k * a);
    Vec3D new_point;
    new_point.x = point.x * cos_z - point.y * sin_z;
    if (fabs(new_point.x) < EPS) {
        new_point.x = 0.;
    }
    new_point.y = point.x * sin_z + point.y * cos_z;
    if (fabs(new_point.y) < EPS) {
        new_point.y = 0.;
    }
    new_point.z = point.z;
    if (fabs(new_point.z) < EPS) {
        new_point.z = 0.;
    }
    return new_point;
}

Vec3D rotateVector(const Vec3D& vec, const Vec3D& rot)
{
    Vec3D new_vec = Rotation_x(vec, rot.x);
    new_vec = Rotation_y(new_vec, rot.y);
    new_vec = Rotation_z(new_vec, rot.z);
    return new_vec;
}

void render(const std::vector<Sphere>& spheres, const Cube& cube, const Vec3D& camera, const Vec3D& rotation, const std::vector<Light>& lights) {
    const int width = 768;
    const int height = 768;
    Pixel* imageBuffer = new Pixel[width * height];

    Vec3D rotated_camera = rotateVector(camera, rotation);

    for (int j = 0; j < height; ++j) {
        for (int i = 0; i < width; ++i) {
            float x = (float)(i - width / 2.) / width;
            float y = (float)(j - height / 2.) / height;
            float z = -48.;
            // float x = -(float)width / 2 + i + 0.5;
            // float y = -(float)height / 2 + j + 0.5;
            // float z = 768;
            Vec3D current_point(x, y, z);
            Vec3D dir = Vec3D((current_point.x - camera.x), (current_point.y - camera.y), (current_point.z - camera.z));
            dir = rotateVector(dir, rotation);
            Ray ray = Ray(rotated_camera, dir);
            imageBuffer[i + j * width] = cast_ray(ray, spheres, cube, lights);
        }
    }
    stbi_flip_vertically_on_write(0);
    stbi_write_png("../out.png", width, height, sizeof(Pixel), imageBuffer, sizeof(Pixel) * width);
}


int main()
{
    std::vector<Sphere> spheres;
    std::vector<Light> lights;

    Material ivory(Pixel(102, 102, 78));
    Material red_rubber(Pixel(178, 26, 26));
    Material dark_orchid(Pixel(153, 50, 204));
    

    Sphere sphere1(Vec3D(5, -6 , 5), 2, ivory);
    Sphere sphere2(Vec3D(0, -6 , 6), 2, dark_orchid);
    Cube cube(Vec3D(0, 0, 0), 4, red_rubber); 
    Vec3D default_camera(0, 0, -50);
    Vec3D rotation(30, 45, 0);
    Light light1(Vec3D(0, -6, 0), 1);
    Light light2(Vec3D(8, -6, 0), 1.2);

    spheres.push_back(sphere1);
    spheres.push_back(sphere2);
    lights.push_back(light1);
    lights.push_back(light2);
    render(spheres, cube, default_camera, rotation, lights);
    return 0;
}