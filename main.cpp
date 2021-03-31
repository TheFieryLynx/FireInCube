#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include "templates.h"
#include "functions.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"


bool triangle_intersection(const Vec3D& orig, const Vec3D& dir, const Vec3D& v0,
                                const Vec3D& v1, const Vec3D& v2) {
    
    // const float EPSILON = 0.0000001;
    // Vec3D edge1, edge2, h, s, q;
    // float a,f,u,v;
    // edge1 = v1 - v0;
    // edge2 = v2 - v0;
    // h = cross(rayVector, edge2);
    // a = dot(edge1, h);
    // if (a > -EPSILON && a < EPSILON)
    //     return false;    // This ray is parallel to this triangle.
    // f = 1.0/a;
    // s = rayOrigin - v0;
    // u = f * dot(s, h);
    // if (u < 0.0 || u > 1.0)
    //     return false;
    // q = cross(s, edge1);
    // v = f * dot(rayVector, q);
    // if (v < 0.0 || u + v > 1.0)
    //     return false;
    // // At this stage we can compute t to find out where the intersection point is on the line.
    // float t = f * dot(edge2,q);
    // if (t > EPSILON) // ray intersection
    // {
    //     Vec3D outIntersectionPoint = rayOrigin + rayVector * t;
    //     return true;
    // }
    // else // This means that there is a line intersection but not a ray intersection.
    //     return false;
    
    Vec3D e1 = v1 - v0;
    Vec3D e2 = v2 - v0;
    // Вычисление вектора нормали к плоскости
    Vec3D pvec = cross(dir, e2);
    float det = dot(e1, pvec);

    // Луч параллелен плоскости
    if (det < 1e-8 && det > -1e-8) {
        return 0;
    }

    float inv_det = 1 / det;
    Vec3D tvec = orig - v0;
    float u = dot(tvec, pvec) * inv_det;
    if (u < 0 || u > 1) {
        return 0;
    }

    Vec3D qvec = cross(tvec, e1);
    float v = dot(dir, qvec) * inv_det;
    if (v < 0 || u + v > 1) {
        return 0;
    }
    return dot(e2, qvec) * inv_det;
}

struct Cube
{
    Vec3D center;
    float radius;
    Rotation rotation;

    float cos_x, sin_x;
    float cos_y, sin_y;
    float cos_z, sin_z;
    Vec3D a, a1;
    Vec3D b, b1;
    Vec3D c, c1;
    Vec3D d, d1;
    Cube(const Vec3D &c_, const float &r_, const Rotation &rot_) : center(c_), radius(r_), rotation(rot_) {
        setCosSin();
        //std::cout << center << std::endl;
        setCoords(a, center.x - radius, center.y + radius, center.z - radius);
        setCoords(b, center.x - radius, center.y + radius, center.z + radius);
        setCoords(c, center.x + radius, center.y + radius, center.z + radius);
        setCoords(d, center.x + radius, center.y + radius, center.z - radius);

        setCoords(a1, center.x - radius, center.y - radius, center.z - radius);
        setCoords(b1, center.x - radius, center.y - radius, center.z + radius);
        setCoords(c1, center.x + radius, center.y - radius, center.z + radius);
        setCoords(d1, center.x + radius, center.y - radius, center.z - radius);
        std::cout << a << std::endl;
        std::cout << a1 << std::endl;
        std::cout << d1 << std::endl;
    }

    void setCoords(Vec3D& a, float x, float y, float z) {
        a.x = x;
        a.y = y;
        a.z = z;
    }

    void setCosSin() {
        float k = 180 / M_PI;
        cos_x = std::cos(k * rotation[0]);
        sin_x = std::sin(k * rotation[0]);
        cos_y = std::cos(k * rotation[1]);
        sin_y = std::sin(k * rotation[1]);
        cos_z = std::cos(k * rotation[2]);
        sin_z = std::sin(k * rotation[2]);
    }

    // Vec3D Rotation_x(Vec3D point) {
    //     Vec3D new_point;
    //     new_point.x = point.x;
    //     new_point.y = point.y + cos_x + point.z * sin_x;
    //     new_point.z = -point.y * sin_x + point.z * cos_x;
    //     return new_point;
    // }

    // Vec3D Rotation_y(Vec3D point) {
    //     Vec3D new_point;
    //     new_point.x = point.x * cos_y + point.z * sin_y;
    //     new_point.y = point.y;
    //     new_point.z = -point.x * sin_y + point.z * cos_y;
    //     return new_point;
    // }

    // Vec3D Rotation_z(Vec3D point) {
    //     Vec3D new_point;
    //     new_point.x = point.x * cos_z - point.y * sin_z;
    //     new_point.y = -point.x * sin_z + point.y * cos_z;
    //     new_point.z = point.z;
    //     return new_point;
    // }

    bool ray_intersect(const Vec3D &orig, const Vec3D &dir) const{
        // std::cout << orig << std::endl;
        // std::cout << dir << std::endl;
        // std::cout << a << std::endl;
        // std::cout << a1 << std::endl;
        // std::cout << d1 << std::endl << std::endl;;
        if (!triangle_intersection(orig, dir, a, a1, d1)) {
            return false;
        }
       
        return true;
    }
};


Pixel cast_ray(const Vec3D &orig, const Vec3D &dir, const Cube &cube) {
    //float cube_dist = std::numeric_limits<float>::max();
    if (!cube.ray_intersect(orig, dir)) {
        //std::cout << "======" << std::endl;
        return Pixel(255 * 0.2, 255 * 0.7, 255 * 0.8); // background color

    }
    std::cout << "." << std::endl;
    return Pixel(255 * 0.4, 255 * 0.4, 255 * 0.3);
}


void render(const Cube &cube) 
{
    const size_t width = 512;
    const size_t height = 512;
    const int fov = M_PI/2.;
    Pixel* buffer = new Pixel[width * height];

    //#pragma omp parallel for
    for (size_t j = 0; j < height; ++j) {
        for(size_t i = 0; i < width; ++i) {
            //buffer[i + j * width] = Pixel(255 * j / float(height), 255 * i / float(width), 0);
            //(float)x*wv/wc/2, (float)y*hv/hc/2, d
            // float x =  (2*(i + 0.5)/(float)width  - 1)*tan(fov/2.)*width/(float)height;
            // float y = -(2*(j + 0.5)/(float)height - 1)*tan(fov/2.);
            float x = - (float)width / 2 + i;
            float y = - (float)height / 2 + j;
            Vec3D dir = Vec3D(x, y, 25).normalize();
            //std::cout << dir << std::endl;
            //std::cout << dir << std::endl;
            buffer[i + j * width] = cast_ray(Vec3D(0, 0, 0), dir, cube);
        }
    }


    stbi_write_png("../out.png", width, height, sizeof(Pixel), buffer, sizeof(Pixel) * width);
}

int main()
{
    Cube cube(Vec3D(0, 0, 50), 100, Rotation(0, 0, 0));
    render(cube);
    std::cout << "Hello World" << std::endl;
    return 0;
}