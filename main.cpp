#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>
#include "geom.h"


#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

int envmap_width, envmap_height;
std::vector<Pixel> envmap;

float cosinus(const Vec3D& a, const Vec3D& b) 
{
    float m = dot(a, b);
    float n1 = std::sqrt(std::pow(a.x, 2) + std::pow(a.y, 2) + std::pow(a.z, 2));
    float n2 = std::sqrt(std::pow(b.x, 2) + std::pow(b.y, 2) + std::pow(b.z, 2));
    return m / (n1 * n2);
}

template <typename Object>
bool scene_intersect(const Ray& ray, Object &obj, float &dist_i, Vec3D& point, Vec3D& N, Material &material) {
    float obj_dist = std::numeric_limits<float>::max();
    if (obj.ray_intersect(ray, dist_i)) {
        obj_dist = dist_i;
        point = ray.orig + ray.dir * dist_i;
        N = obj.normInPoint(point);
        material = obj.material;
        //std::cout << material.diffuse_color << std::endl;
    }
    return obj_dist < 10000;
    // float checkerboard_dist = std::numeric_limits<float>::max();
    // if (fabs(ray.dir.y)>1e-3)  {
    //     float d = -(ray.orig.y+4)/ray.dir.y; // the checkerboard plane has equation y = -4
    //     Vec3D pt = ray.orig + ray.dir*d;
    //     if (d > 0 && fabs(pt.x)<10 && pt.z<-10 && pt.z>-30 && d<spheres_dist) {
    //         checkerboard_dist = d;
    //         point = pt;
    //         N = Vec3D(0,1,0);
    //         material.diffuse_color = (int(.5*point.x+1000) + int(.5*point.z)) & 1 ? Pixel(1,1,1) : Pixel(1, .7, .3);
    //         material.diffuse_color = material.diffuse_color * .3;
    //     }
    // }
    //std::cout << checkerboard_dist << std::endl;
    //return std::min(spheres_dist, checkerboard_dist)<1000;
}

template <typename Object>
bool cast_ray_object(const Ray& ray, WorldObjects& world_obj, Object &obj, float &dist_i, Vec3D& point, Vec3D& N, Material &material)
{
    if(scene_intersect(ray, obj, dist_i, point, N, material)) {
        //std::cout << "!" << std::endl;
        float diffuse_light_intensity = 0.3;
        float specular_light_intensity = 0;
        for (int i = 0; i < int(world_obj.lights.size()); ++i) {
            Vec3D light_dir = (world_obj.lights[i].position - point).normalize();
            bool check = true;

            float light_distance = (world_obj.lights[i].position - point).norm();
            //std::cout << N << std::endl;
            Vec3D shadow_orig = dot(light_dir, N) < 0 ? point - N * 0.001 : point + N * 0.001; // checking if the point lies in the shadow of the world_obj.lights[i]
            //std::cout << "#" << std::endl;
            Vec3D tmp_pt, shadow_pt, shadow_N;
            float t_tmp, t = std::numeric_limits<float>::max();
            Ray new_ray(shadow_orig, light_dir);
            Material tmpmaterial;
            for(auto j : world_obj.spheres) {
                if (scene_intersect(new_ray, j, t_tmp, tmp_pt, shadow_N, tmpmaterial)) {
                    if (t_tmp < t) {
                        shadow_pt = tmp_pt;
                        t = t_tmp;
                    }
                    if ((shadow_pt - shadow_orig).norm() < light_distance) {
                        check = false;
                    }
                }
            }
            if (check) {
                if (scene_intersect(new_ray, world_obj.cubes[0], t_tmp, tmp_pt, shadow_N, tmpmaterial)) {
                    if (t_tmp < t) {
                        shadow_pt = tmp_pt;
                        t = t_tmp;
                    }
                    if ((shadow_pt - shadow_orig).norm() < light_distance) {
                        check = false;
                    }
                }
            }

            if (check) {
                diffuse_light_intensity  += world_obj.lights[i].intensity * std::max(0.f, dot(N, light_dir));  //fabs(cosinus(light_dir, N)); //std::max(0.f, light_dir * N);
                specular_light_intensity += powf(std::max(0.f, dot(reflect(light_dir, N), -ray.dir)), material.specular_exponent) * world_obj.lights[i].intensity;
            }  
        }
        Pixel pix(1., 1., 1.);
        material.diffuse_color = material.diffuse_color * diffuse_light_intensity * material.albedo.x + specular_light_intensity * material.albedo.y;;
        return true;
    }
    return false;
}

/*
bool cast_ray_sphere(const Ray& ray, const WorldObjects& world_obj, const Sphere &sphere, float &dist_i, Vec3D& point, Vec3D& N, Material &material) 
{   
    if(scene_intersect(ray, sphere, dist_i, point, N, material)) {
        float diffuse_light_intensity = 0.5;
        float specular_light_intensity = 0;
        
        for (int i = 0; i < int(world_obj.lights.size()); ++i) {
            Vec3D light_dir = (world_obj.lights[i].position - point).normalize();
            bool check = true;
            
            // float light_distance = (world_obj.lights[i].position - point).norm();
            // //std::cout << N << std::endl;
            // Vec3D shadow_orig = dot(light_dir, N) < 0 ? point - N * 0.001 : point + N * 0.001; // checking if the point lies in the shadow of the world_obj.lights[i]
            // //std::cout << "#" << std::endl;
            // Vec3D tmp_pt, shadow_pt, shadow_N;
            // float t_tmp, t = std::numeric_limits<float>::max();
            // Ray new_ray(shadow_orig, light_dir);
            // Material tmpmaterial;
            // for(auto j : world_obj.spheres) {
            //     if (scene_intersect(new_ray, j, t_tmp, tmp_pt, shadow_N, tmpmaterial)) {
            //         //check = false;
            //         if (t_tmp < t) {
            //             shadow_pt = tmp_pt;
            //             t = t_tmp;
            //         }
            //         if ((shadow_pt - shadow_orig).norm() < light_distance) {
            //             //std::cout << shadow_pt << " - " << shadow_orig << " < " << light_distance << std::endl;
            //             check = false;
            //         }
            //     }
            // }

            if (check) {
                diffuse_light_intensity  += world_obj.lights[i].intensity * std::max(0.f, dot(N, light_dir));  //fabs(cosinus(light_dir, N)); //std::max(0.f, light_dir * N);
                specular_light_intensity += powf(std::max(0.f, dot(reflect(light_dir, N), -ray.dir)), material.specular_exponent) * world_obj.lights[i].intensity;
            }   
            
        }
        
        material.diffuse_color = material.diffuse_color * diffuse_light_intensity * material.albedo.x + Pixel(1., 1., 1.) * specular_light_intensity * material.albedo.y;;
        return true;
    }
    return false;
}

bool cast_ray_cube(const Ray& ray, const Cube &cube, float& t, Vec3D& point, Vec3D& N, Material &material, const std::vector<Light>& lights) 
{
    if(scene_intersect(ray, cube, t, point, N, material)) {
        float diffuse_light_intensity = 0.5;
        float specular_light_intensity = 0;
        for (int i = 0; i < int(lights.size()); ++i) {
            Vec3D light_dir = (lights[i].position - point).normalize();
            diffuse_light_intensity  += lights[i].intensity * std::max(0.f, dot(N, light_dir));  //fabs(cosinus(light_dir, N)); //std::max(0.f, light_dir * N);
            specular_light_intensity += powf(std::max(0.f, dot(reflect(light_dir, N), -ray.dir)), material.specular_exponent) * lights[i].intensity;
        }
        material.diffuse_color = material.diffuse_color * diffuse_light_intensity * material.albedo.x + Pixel(1., 1., 1.) * specular_light_intensity * material.albedo.y;;
        return true;
    }
    return false;
}
*/

Pixel cast_ray(const Ray& ray, WorldObjects& world_obj, size_t depth)
{
    Vec3D point, N;
    std::vector<int> flags;
    for (int i = 0; i < int(world_obj.spheres.size()) + 2; ++i) flags.push_back(1);
    std::vector<float> t_s;
    std::vector<Material> materials;

    float t_cube = std::numeric_limits<float>::max();
    Material cube_material;
    if (depth > 4 || !cast_ray_object(ray, world_obj, world_obj.cubes[0], t_cube, point, N, cube_material)) {
        flags[0] = 0;
        t_cube = -1.0;
    }
    materials.push_back(cube_material);
    t_s.push_back(t_cube);

    float t_surface = std::numeric_limits<float>::max();
    Material surface_material;
    if (depth > 4 || !cast_ray_object(ray, world_obj, world_obj.surfaces[0], t_surface, point, N, surface_material)) {
        //std::cout << "@" << std::endl;
        flags[1] = 0;
        t_surface = -1.0;
    }
    //flags[1] = 0;
    //std::cout << surface_material.diffuse_color << std::endl;
    materials.push_back(surface_material);
    t_s.push_back(t_surface);

    for (int num_s = 0; num_s < int(world_obj.spheres.size()); ++num_s) {
        float t_sphere = std::numeric_limits<float>::max();
        Material sphere_material;
        if (depth > 4 || !cast_ray_object(ray, world_obj, world_obj.spheres[num_s], t_sphere, point, N, sphere_material)) {
            flags[num_s + 2] = 0;
            t_sphere = -1.0;
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

    Vec3f reflect_dir = reflect(dir, N).normalize();

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

void render(WorldObjects& world_obj, const Vec3D& camera, const Vec3D& rotation) {
    const int width = 1024;
    const int height = 1024;
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
            imageBuffer[i + j * width] = cast_ray(ray, world_obj, 0);
        }
    }
    //std::vector<unsigned char> pixmap(width*height*3);

    for (size_t i = 0; i < height * width; ++i) {
        //Pixel &c = imageBuffer[i];
        float max = std::max(imageBuffer[i].R, std::max(imageBuffer[i].G, imageBuffer[i].B));
        if (max > 1.) {
            imageBuffer[i] = imageBuffer[i] / max;
        }
        //Гамма коорекция
        imageBuffer[i].R = std::powf(imageBuffer[i].R, 1. / 2.2);
        imageBuffer[i].G = std::powf(imageBuffer[i].G, 1. / 2.2);
        imageBuffer[i].B = std::powf(imageBuffer[i].B, 1. / 2.2);


        // pixmap[i*3 + 0] = (unsigned char)(255 * std::max(0.f, std::min(1.f, imageBuffer[i].R)));
        // pixmap[i*3 + 1] = (unsigned char)(255 * std::max(0.f, std::min(1.f, imageBuffer[i].G)));
        // pixmap[i*3 + 2] = (unsigned char)(255 * std::max(0.f, std::min(1.f, imageBuffer[i].B)));
        finalBuffer[i] = Pixel_final(imageBuffer[i].R * 255, imageBuffer[i].G * 255, imageBuffer[i].B * 255);
    }
    stbi_flip_vertically_on_write(0);
    stbi_write_png("../out.png", width, height, sizeof(Pixel_final), finalBuffer, sizeof(Pixel_final) * width);
}

int main()
{
    WorldObjects world_obj;

    Material ivory(Vec3D(0.6, 0.3, 0.), Pixel(0.4, 0.4, 0.3), 50.);
    Material red_rubber(Vec3D(0.9, 0.1, 0.), Pixel(0.3, 0.1, 0.1), 50.);
    Material dark_orchid(Vec3D(0.9, 0.5, 0.), Pixel(0.6, 0.19, 0.8), 50.);

    Material dark_gray(Vec3D(0.9, 0.1, 0.), Pixel(0.1, 0.4, 0.22), 50.);
    Material light_gray(Vec3D(0.9, 0.1, 0.), Pixel(0.15, 1., 0.34), 50.);
    
    Material mirror(Vec3D(0.0, 10.0, 0.8), Pixel(1.0, 1.0, 1.0), 1425.);

    Surface surface(4, 1.5, dark_gray, light_gray);

    Sphere sphere1(Vec3D(5, -6 , 5), 2, ivory);
    Sphere sphere2(Vec3D(0, -6 , 6), 2, red_rubber);
    Cube cube(Vec3D(0, 0, 0), 4, ivory); 

    Vec3D default_camera(0, 0, -50);
    Vec3D rotation(10, 35, 0);

    Light light1(Vec3D(0, -6, 0), 1.2);
    Light light2(Vec3D(10, -6, 1), 1);
    Light light3(Vec3D(50, -20, -50), 1.2);
    Light light4(Vec3D(100, -15 , 100), 1.8);
    Light light5(Vec3D(5, -2, -6), 1);

    world_obj.surfaces.push_back(surface);

    world_obj.spheres.push_back(sphere1);
    world_obj.spheres.push_back(sphere2);

    world_obj.cubes.push_back(cube);

    world_obj.lights.push_back(light1);
    world_obj.lights.push_back(light2);
    world_obj.lights.push_back(light3);
    world_obj.lights.push_back(light4);
    //world_obj.lights.push_back(light5);

    render(world_obj, default_camera, rotation);
    return 0;
}