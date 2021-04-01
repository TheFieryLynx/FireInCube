#ifndef __TEMPLATES_H__
#define __TEMPLATES_H__

struct Vec3D
{
    Vec3D() {
        x = 0.;
        y = 0.;
        z = 0.;
    }
    Vec3D(const float& x_, const float& y_, const float& z_) : x(x_), y(y_), z(z_) {}
    float norm();
    Vec3D normalize();
    float x;
    float y;
    float z;
};

Vec3D operator-(const Vec3D& l_v, const Vec3D& r_v)
{
    float x = l_v.x - r_v.x;
    float y = l_v.y - r_v.y;
    float z = l_v.z - r_v.z;
    return Vec3D(x, y, z);
}

Vec3D operator*(const Vec3D& l_v, float c) 
{
    float x = l_v.x * c;
    float y = l_v.y * c;
    float z = l_v.z * c;
    return Vec3D(x, y , z);
}

float Vec3D::norm() {
    return std::sqrt(x * x + y * y + z * z);
}

Vec3D Vec3D::normalize() {
    *this = (*this) * (1. / norm());
    return *this;
}

struct Pixel
{
    Pixel(){}
    Pixel(const int& r_, const int& g_, const int& b_) : R(r_), G(g_), B(b_) {}
    uint8_t R;
    uint8_t G;
    uint8_t B;
};

struct Light {
    Light(const Vec3D& p, const float& i) : position(p), intensity(i) {}
    Vec3D position;
    float intensity;
 };

class Ray 
{ 
public: 
    Ray(const Vec3D &orig, const Vec3D &dir) : orig(orig), dir(dir) 
    { 
        invdir = Vec3D(1 / dir.x, 1 / dir.y, 1 / dir.z);
        sign[0] = (invdir.x < 0); 
        sign[1] = (invdir.y < 0); 
        sign[2] = (invdir.z < 0); 
    } 
    Vec3D orig, dir; // ray orig and dir 
    Vec3D invdir; 
    int sign[3]; 
}; 
 
class Cube 
{ 
public: 
    Cube(const Vec3D &b0, const Vec3D &b1) { 
        bounds[0] = b0; 
        bounds[1] = b1; 
    } 
    bool ray_intersect(const Ray&, float&) const;
    Vec3D bounds[2]; 
}; 

bool Cube::ray_intersect(const Ray &r, float &t) const 
{ 
    float tmin, tmax, tymin, tymax, tzmin, tzmax; 

    tmin = (bounds[r.sign[0]].x - r.orig.x) * r.invdir.x; 
    tmax = (bounds[1-r.sign[0]].x - r.orig.x) * r.invdir.x; 
    tymin = (bounds[r.sign[1]].y - r.orig.y) * r.invdir.y; 
    tymax = (bounds[1-r.sign[1]].y - r.orig.y) * r.invdir.y; 

    if ((tmin > tymax) || (tymin > tmax)) 
        return false; 

    if (tymin > tmin) 
    tmin = tymin; 
    if (tymax < tmax) 
    tmax = tymax; 

    tzmin = (bounds[r.sign[2]].z - r.orig.z) * r.invdir.z; 
    tzmax = (bounds[1-r.sign[2]].z - r.orig.z) * r.invdir.z; 

    if ((tmin > tzmax) || (tzmin > tmax)) 
        return false; 

    if (tzmin > tmin) 
    tmin = tzmin; 
    if (tzmax < tmax) 
    tmax = tzmax; 

    t = tmin; 

    if (t < 0) { 
        t = tmax; 
        if (t < 0) return false; 
    } 
    return true; 
} 

Vec3D cross (const Vec3D& a, const Vec3D& b) 
{
    float x = a.y * b.z - a.z * b.y;
    float y = a.z * b.x - a.x * b.z;
    float z = a.x * b.y - a.y * b.x;
    return Vec3D(x, y, z);
}

float dot(const Vec3D& a, const Vec3D &b) 
{   
    return a.x * b.x + a.y * b.y + a.z + b.z;
}

#endif //__TEMPLATES_H__