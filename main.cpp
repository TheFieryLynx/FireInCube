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

/* Vec3D refract(const Vec3D &I, Vec3D N, float outcoming,  float incoming = 1.f) {
    // if (dot(I, N) < 0) {
    //     N = -N;
    // } else {
    //     float tmp = outcoming;
    //     outcoming = incoming;
    //     incoming = tmp;
    // }
    // float cos_t1 = cosinus(I, N);
    // std::cout << "I: " << I << std::endl;
    // std::cout << "N: " << N << std::endl;
    // std::cout << "cos: " << cos_t1 << std::endl;
    // float sin_t1 = std::sqrt(1.0 - cos_t1 * cos_t1);
    // Vec3D tang = (I - N * cos_t1).normalize();
    // std::cout << incoming << " " << outcoming << std::endl;
    // float n1n2 = incoming / outcoming;
    // float sin_t2 = sin_t1 * n1n2;
    // if (sin_t2 > 1) {
    //     return reflect(-I, N);
    // }
    // float cos_t2 = std::sqrt (1 - sin_t2 * sin_t2);
    // return tang * sin_t2 + N * cos_t2;

    //std::cout << "I: " << I << std::endl;
    //std::cout << "N: " << N << std::endl;
    if (cosinus(I, N) < 0) {
        //std::cout << "!" << std::endl;
        N = -N;
        float tmp = outcoming;
        outcoming = incoming;
        incoming = tmp;
    } else {
        //std::cout << "!!" << std::endl;
        float tmp = outcoming;
        outcoming = incoming;
        incoming = tmp;
    }
    float cos_A = cosinus(I, N);
    //std::cout << "cos_A: " << cos_A << std::endl;

    float n1n2 = incoming / outcoming;
    //std::cout << incoming << " " << outcoming << std::endl;

    float sin_B = n1n2 * std::sqrtf(1 - cos_A * cos_A);
    //std::cout << "sin_B: " << sin_B << std::endl;
    if (sin_B > 1) {
        return reflect(-I, N);
    }  
    float cos_B = std::sqrtf(1 - sin_B * sin_B);
    return I * n1n2 + N * (n1n2 * cos_A - cos_B);
} */

Vec3D refract(const Vec3D &I, const Vec3D &N, const float outcoming, const float incoming=1.f) { 
    float cosi = -cosinus(I, N);
    if (cosi < 0) {
        return refract(I, -N, incoming, outcoming);
    } 
       
    float eta = incoming / outcoming;

    float k = 1 - eta*eta*(1 - cosi*cosi);

    return k < 0 ? reflect(-I, N) : I*eta + N*(eta*cosi - sqrtf(k)); 
}

template <typename Object>
void object_intersect(const Ray& ray, std::vector<Object>& objs, float &obj_dist, Vec3D& point, Vec3D& N, Material &material) {
    float dist_i;
    for(int i = 0; i < int(objs.size()); ++i) {
        if (objs[i].ray_intersect(ray, dist_i) && dist_i < obj_dist) {
            obj_dist = dist_i;
            point = ray.orig + ray.dir * dist_i;
            N = objs[i].normInPoint(point);
            material = objs[i].material;
            //std::cout << material.diffuse_color.G << std::endl;
        }
    }
}

bool scene_intersect(const Ray& ray, WorldObjects& world_obj, Vec3D& point, Vec3D& N, Material &material) {
    float obj_dist = std::numeric_limits<float>::max();
    
    object_intersect(ray, world_obj.spheres, obj_dist, point, N, material);
    object_intersect(ray, world_obj.surfaces, obj_dist, point, N, material);
   
    object_intersect(ray, world_obj.cubes, obj_dist, point, N, material);
    object_intersect(ray, world_obj.cones, obj_dist, point, N, material);
    return obj_dist < 10000;
}



Pixel cast_ray(const Ray& ray, WorldObjects& world_obj, size_t depth)
{   
    static float current_reflactive_index = 1.0;
    static float previous_reflactive_index = 1.0;
    Vec3D point, N;
    Material material;
    cube_intersect_bottom = false;
    if(depth < 5 && scene_intersect(ray, world_obj, point, N, material)) {
        
        Vec3D reflect_dir = reflect(-ray.dir, N).normalize();
        //std::cout << material.refractive_index << std::endl;
        /*
            if (material.refractive_index != refract_ind[i]) {
            refract_ind[i + 1] = material.refractive_index;
            current_reflactive_index = refract_ind[i + 1];
            previous_reflactive_index = refract_ind[i];
            ++i;
            if (i == 3) {
                i = 2;
            }
        } else {
            --i;
            current_reflactive_index = refract_ind[i + 1];
            previous_reflactive_index = refract_ind[i];
            refract_ind[i + 1] = 1.0;
        }
        */
        
        Vec3D refract_dir = refract(ray.dir, N, material.refractive_index).normalize();

        Vec3D reflect_orig = dot(reflect_dir, N) < 0 ? point - N * 1e-3 : point + N * 1e-3;
        Vec3D refract_orig = dot(refract_dir, N) < 0 ? point - N * 1e-3 : point + N * 1e-3;


        Pixel reflect_color = cast_ray(Ray(reflect_orig, reflect_dir), world_obj, depth + 1);
        Pixel refract_color = cast_ray(Ray(refract_orig, refract_dir), world_obj, depth + 1);
        
        
        float diffuse_light_intensity = 0.2;
        float specular_light_intensity = 0;
        for (int i = 0; i < int(world_obj.lights.size()); ++i) {
            Vec3D light_dir = (world_obj.lights[i].position - point).normalize();
            bool check = true;
            
            float light_distance = (world_obj.lights[i].position - point).norm();
            //std::cout << N << std::endl;
            Vec3D shadow_orig = dot(light_dir, N) < 0 ? point - N * 0.001 : point + N * 0.001; // checking if the point lies in the shadow of the world_obj.lights[i]
            //std::cout << "#" << std::endl;
            Vec3D tmp_pt, shadow_pt, shadow_N;
            //float t = std::numeric_limits<float>::max();
            Ray new_ray(shadow_orig, light_dir);
            Material material_tmp;

            if (scene_intersect(new_ray, world_obj, shadow_pt, shadow_N, material_tmp)) {
                if ((shadow_pt - shadow_orig).norm() < light_distance) {
                    check = false;
                }
            }
            
            if (check) {
                diffuse_light_intensity  += world_obj.lights[i].intensity * std::max(0.f, dot(N, light_dir));
                specular_light_intensity += powf(std::max(0.f, dot(reflect(light_dir, N), -ray.dir)), material.specular_exponent) * world_obj.lights[i].intensity;
            }  
        }
        Pixel pix(1., 1., 1.);
        material.diffuse_color = material.diffuse_color * diffuse_light_intensity * material.albedo.x + 
                            pix * specular_light_intensity * material.albedo.y + reflect_color * material.albedo.z + refract_color * material.albedo.a;

        return material.diffuse_color;
    }
    return Pixel(0., 0.78, 1.);
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

/*Pixel cast_ray(const Ray& ray, WorldObjects& world_obj)
{
    Vec3D point, N;
    std::vector<int> flags;
    for (int i = 0; i < int(world_obj.spheres.size()) + 1; ++i) flags.push_back(1);
    std::vector<float> t_s;
    std::vector<Pixel> colors;

    float t_cube = std::numeric_limits<float>::max();
    Material cube_material;
    Pixel pix = cast_ray_object(ray, world_obj, t_cube, 0);
    if (pix.R == 0 && pix.G == 0.78 && pix.R == 1.) {
        flags[0] = 0;
        t_cube = -1.0;
    }
    colors.push_back(pix);
    t_s.push_back(t_cube);

    // float t_surface = std::numeric_limits<float>::max();
    // Material surface_material;
    // pix = cast_ray_object(ray, world_obj, world_obj.surfaces[0], t_surface, 0);
    // if (pix.R == 0 && pix.G == 0.78 && pix.R == 1.) {
    //     //std::cout << "@" << std::endl;
    //     flags[1] = 0;
    //     t_surface = -1.0;
    // }
    // //flags[1] = 0;
    // //std::cout << surface_material.diffuse_color << std::endl;
    // colors.push_back(pix);
    // t_s.push_back(t_surface);

    for (int num_s = 0; num_s < int(world_obj.spheres.size()); ++num_s) {
        float t_sphere = std::numeric_limits<float>::max();
        Material sphere_material;
        pix = cast_ray_object(ray, world_obj, t_sphere, 0);
        if (pix.R == 0 && pix.G == 0.78 && pix.R == 1.) {
            flags[num_s + 1] = 0;
            t_sphere = -1.0;
        }
        colors.push_back(pix);
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
    return colors[i_tmp];
} */

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

    Material ivory(1.0, Vec4D(0.6, 0.3, 0., 0.), Pixel(0.6, 0.4, 0.3), 50.);
    Material red_rubber(1.0, Vec4D(0.9, 0.1, 0., 0.), Pixel(0.3, 0.5, 0.3), 50.);
    Material dark_orchid(1.0, Vec4D(0.9, 0.5, 0., 0.), Pixel(0.53, 0.1, 0.9), 50.);

    Material dark_gray(1.0, Vec4D(0.9, 0.1, 0., 0.), Pixel(0.44, 0.44, 0.44), 50.);
    Material light_gray(1.0, Vec4D(0.9, 0.1, 0., 0.), Pixel(0.22, 0.22, 0.22), 50.);
    
    Material mirror(1.0, Vec4D(0.0, 10.0, 0.8, 0.), Pixel(1.0, 1.0, 1.0), 100.);

    Material glass(1.3, Vec4D(0.0, 0.5, 0.1, 0.8), Pixel(0.6, 0.7, 0.8), 125.);

    Surface surface(5, 1.5, dark_gray, light_gray);
    //---------------------------------<Objects>---------------------------------;
    Sphere sphere1(Vec3D(5, -6 , 5), 2, ivory);
    Sphere sphere2(Vec3D(0, -6 , 6), 2, red_rubber);

    Sphere sphere3(Vec3D(-8, 0 , 1), 2, glass);
    Sphere sphere4(Vec3D(0, 0, 10), 2, red_rubber);


    //Sphere sphere_in_cube(Vec3D(0, 0, 0), 2, ivory);

    Cube cube(Vec3D(0, 0, 0), 4, glass); 

    Cone cone(Vec3D(0, 5., 0), 8, 1, ivory, dark_orchid);

    //---------------------------------<Camera>---------------------------------;

    Vec3D default_camera(0, 0, -50);
    Vec3D rotation(20, 30, 0);

    //---------------------------------<Lights>---------------------------------;
    Light light1(Vec3D(0, -6, 0), 1.2);
    Light light2(Vec3D(10, -6, 1), 1);
    Light light3(Vec3D(50, -20, -50), 1.2);
    Light light4(Vec3D(100, -15 , 100), 1.2);
    Light light5(Vec3D(2, -8 , -20), 1.5);

    world_obj.surfaces.push_back(surface);

    //world_obj.spheres.push_back(sphere1);
    //world_obj.spheres.push_back(sphere2);
    world_obj.spheres.push_back(sphere3);
    world_obj.spheres.push_back(sphere4);
    //world_obj.spheres.push_back(sphere5);

    world_obj.cubes.push_back(cube);

    world_obj.cones.push_back(cone);

    world_obj.lights.push_back(light1);
    world_obj.lights.push_back(light2);
    // world_obj.lights.push_back(light3);
    // world_obj.lights.push_back(light4);
    //world_obj.lights.push_back(light5);
    std::cout << "----<start>----" << std::endl;
    render(world_obj, default_camera, rotation);
    Vec3D newdir = refract(Vec3D(1, 0, 1), Vec3D(-1, 0, 0), 1.);
    std::cout << "Inside: " <<newdir << std::endl << std::endl;
    Vec3D newdir2 = refract(newdir, Vec3D(0, 0, 1), 1.);
    std::cout << "Outside: " << newdir2 << std::endl;
    return 0;
}