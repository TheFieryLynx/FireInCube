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

bool scene_sphere_intersect(const Ray& ray, const Sphere &sphere, float &dist_i, Vec3D& point, Vec3D& N, Material &material) {
    if (sphere.ray_intersect(ray.orig, ray.dir, dist_i)) {
        point = ray.orig + ray.dir * dist_i;
        N = (point - sphere.center).normalize();
        material = sphere.material;
        return true;
    }
    return false;
}

bool cast_ray_sphere(const Ray& ray, const Sphere &sphere, float &dist_i, Vec3D& point, Vec3D& N, Material &material, const std::vector<Light>& lights) 
{   
    if(scene_sphere_intersect(ray, sphere, dist_i, point, N, material)) {
        float diffuse_light_intensity = 0;
        float specular_light_intensity = 0;
        
        for (int i = 0; i < int(lights.size()); ++i) {
            Vec3D light_dir = (lights[i].position - point).normalize();

            // float light_distance = (lights[i].position - point).norm();
            // std::cout << "!" << std::endl;
            // Vec3D shadow_orig = dot(light_dir, N) < 0 ? point - N * 1e-3 : point + N * 1e-3; // checking if the point lies in the shadow of the lights[i]
            // std::cout << "#" << std::endl;
            // Vec3D shadow_pt, shadow_N;
            // float t_tmp;
            // Ray new_ray(shadow_orig, light_dir);
            // Material tmpmaterial;
            // if (scene_sphere_intersect(new_ray, sphere, t_tmp, shadow_pt, shadow_N, tmpmaterial, lights) && (shadow_pt - shadow_orig).norm() < light_distance) {
            //     continue;
            // }
            //std::cout << material.specular_exponent << std::endl;
            diffuse_light_intensity  += lights[i].intensity * std::max(0.f, dot(N, light_dir));  //fabs(cosinus(light_dir, N)); //std::max(0.f, light_dir * N);
            specular_light_intensity += powf(std::max(0.f, dot(reflect(light_dir, N), -ray.dir)), material.specular_exponent) * lights[i].intensity;
        }
        
        material.diffuse_color = material.diffuse_color * diffuse_light_intensity * material.albedo.x + Pixel(1., 1., 1.) * specular_light_intensity * material.albedo.y;;
        return true;
    }
    return false;
}

bool scene_cube_intersect(const Ray& ray, const Cube &cube, float& t, Vec3D& point, Vec3D& N, Material &material) 
{
    if(cube.ray_intersect(ray, t)) {
        point = ray.orig + ray.dir * t;
        N = cube.normInPoint(point);
        //std::cout << N << std::endl;
        material = cube.material;
        return true;
    }
    return false;
}

bool cast_ray_cube(const Ray& ray, const Cube &cube, float& t, Vec3D& point, Vec3D& N, Material &material, const std::vector<Light>& lights) 
{
    if(scene_cube_intersect(ray, cube, t, point, N, material)) {
        float diffuse_light_intensity = 0;
        float specular_light_intensity = 0;
        
        for (int i = 0; i < int(lights.size()); ++i) {
            Vec3D light_dir = (lights[i].position - point).normalize();
            diffuse_light_intensity  += lights[i].intensity * std::max(0.f, dot(N, light_dir));  //fabs(cosinus(light_dir, N)); //std::max(0.f, light_dir * N);
            specular_light_intensity += powf(std::max(0.f, cosinus(reflect(light_dir, N), -ray.dir)), material.specular_exponent) * lights[i].intensity;
        }
        material.diffuse_color = material.diffuse_color * diffuse_light_intensity * material.albedo.x + Pixel(1., 1., 1.) * specular_light_intensity * material.albedo.y;;
        return true;
    }
    return false;
}

Pixel cast_ray(const Ray& ray, const std::vector<Sphere>& spheres, const Cube& cube, const std::vector<Light>& lights)
{
    Vec3D point, N;
    std::vector<int> flags;
    for (int i = 0; i < int(spheres.size()) + 1; ++i) flags.push_back(1);
    std::vector<float> t_s;
    std::vector<Material> materials;

    float t_cube = std::numeric_limits<float>::max();
    Material cube_material;
    if (!cast_ray_cube(ray, cube, t_cube, point, N, cube_material, lights)) {
        flags[0] = 0;
        t_cube = -1.0;
        //cube_material = Material(Vec2D(1, 0), Pixel(0., 0.78, 1.), 10.);
        //return Pixel(0, 200, 255); // background color
    }
    materials.push_back(cube_material);
    t_s.push_back(t_cube);

    for (int num_s = 0; num_s < int(spheres.size()); ++num_s) {
        float t_sphere = std::numeric_limits<float>::max();
        Material sphere_material;
        if (!cast_ray_sphere(ray, spheres[num_s], t_sphere, point, N, sphere_material, lights)) {
            flags[num_s + 1] = 0;
            t_sphere = -1.0;
            //sphere_material = Material(Vec2D(1, 0), Pixel(0., 0.78, 1.), 10.);
            //return Pixel(0, 200, 255); // background color
        }
        materials.push_back(sphere_material);
        t_s.push_back(t_sphere);
    }
    float tmp = std::numeric_limits<float>::max();
    int i_tmp = -1;
    for (int i = 0; i < int(flags.size()); ++i) {
        if (flags[i] != 0) {
            if (t_s[i] < tmp) {
                tmp = t_s[i];
                i_tmp = i;
            }
        }
    }
    if (i_tmp == -1) {
        return Pixel(0., 0.78, 1.);
    }
    return materials[i_tmp].diffuse_color;
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
    Pixel_final *finalBuffer = new Pixel_final[width * height];
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
            Vec3D dir = Vec3D((current_point.x - camera.x), (current_point.y - camera.y), (current_point.z - camera.z)).normalize();
            dir = rotateVector(dir, rotation);
            Ray ray = Ray(rotated_camera, dir);
            imageBuffer[i + j * width] = cast_ray(ray, spheres, cube, lights);
        }
    }

    for (size_t i = 0; i < height * width; ++i) {
        Pixel &c = imageBuffer[i];
        float max = std::max(c.R, std::max(c.G, c.B));
        if (max > 1.) {
            c = c * (1. / max);
        }
        finalBuffer[i] = Pixel_final(c.R * 255, c.G * 255, c.B * 255);
    }
    stbi_flip_vertically_on_write(0);
    stbi_write_png("../out.png", width, height, sizeof(Pixel_final), finalBuffer, sizeof(Pixel_final) * width);

    std::ofstream ofs; // save the framebuffer to file
    ofs.open("../out.ppm");
    ofs << "P6\n" << width << " " << height << "\n255\n";
    for (size_t i = 0; i < height*width; ++i) {
        Pixel &c = imageBuffer[i];
        float max = std::max(c.R, std::max(c.G, c.B));
        if (max>1) c = c /max;
        ofs << (char)(255 * std::max(0.f, std::min(1.f, imageBuffer[i].R)));
        ofs << (char)(255 * std::max(0.f, std::min(1.f, imageBuffer[i].G)));
        ofs << (char)(255 * std::max(0.f, std::min(1.f, imageBuffer[i].B)));
    }
    ofs.close();
}


int main()
{
    std::vector<Sphere> spheres;
    std::vector<Light> lights;

    Material ivory(Vec2D(0.6, 0.3), Pixel(0.4, 0.4, 0.3), 50.);
    Material red_rubber(Vec2D(0.9, 0.1), Pixel(0.3, 0.1, 0.1), 50.);
    Material dark_orchid(Vec2D(0.9, 0.5), Pixel(0.6, 0.19, 0.8), 50.);
    

    Sphere sphere1(Vec3D(5, -6 , 5), 2, ivory);
    Sphere sphere2(Vec3D(0, -6 , 6), 2, red_rubber);
    Cube cube(Vec3D(0, 0, 0), 4, ivory); 

    Vec3D default_camera(0, 0, -50);
    Vec3D rotation(30, 45, 0);

    Light light1(Vec3D(0, -5, 0), 1.5);
    Light light2(Vec3D(8, -6, 0), 1.);
    Light light3(Vec3D(-25, -20, -50), 1.4);

    spheres.push_back(sphere1);
    spheres.push_back(sphere2);

    lights.push_back(light1);
    lights.push_back(light2);
    lights.push_back(light3);
    render(spheres, cube, default_camera, rotation, lights);
    return 0;
}