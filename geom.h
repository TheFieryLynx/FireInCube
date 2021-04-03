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

struct Pixel
{
    Pixel(){}
    Pixel(const uint8_t& r_, const uint8_t& g_, const uint8_t& b_) : R(r_), G(g_), B(b_) {}
    uint8_t R;
    uint8_t G;
    uint8_t B;
};

struct Material {
    Material(const Pixel& color) : diffuse_color(color) {}
    Material() : diffuse_color() {}
    Pixel diffuse_color;
};

Vec3D cross (const Vec3D& a, const Vec3D& b); 

Vec3D operator-(const Vec3D& l_v, const Vec3D& r_v)
{
    float x = l_v.x - r_v.x;
    float y = l_v.y - r_v.y;
    float z = l_v.z - r_v.z;
    return Vec3D(x, y, z);
}

Vec3D operator+(const Vec3D& l_v, const Vec3D& r_v)
{
    float x = l_v.x + r_v.x;
    float y = l_v.y + r_v.y;
    float z = l_v.z + r_v.z;
    return Vec3D(x, y, z);
}

Vec3D operator*(const Vec3D& l_v, float c) 
{
    float x = l_v.x * c;
    float y = l_v.y * c;
    float z = l_v.z * c;
    return Vec3D(x, y , z);
}

float operator*(const Vec3D& l_v, const Vec3D& r_v) {
    float ret = l_v.x * r_v.x + l_v.y * r_v.y + l_v.z * r_v.z;
    return ret;
}

float Vec3D::norm() {
    return std::sqrt(x * x + y * y + z * z);
}

Vec3D Vec3D::normalize() {
    *this = (*this) * (1. / norm());
    return *this;
}

std::ostream& operator<<(std::ostream& out, const Vec3D& v) {
    out << v.x << " " << v.y << " " << v.z;;
    return out ;
}

Pixel operator*(const Pixel& l_v, float c) 
{
    uint8_t r, g, b;
    if (l_v.R * c > 255) {
        r = 255;
    } else {
        r = uint8_t(l_v.R * c);
    }

    if (l_v.G * c > 255) {
        g = 255;
    } else {
        g = uint8_t(l_v.G * c);
    }

    if (l_v.B * c > 255) {
        b = 255;
    } else {
        b = uint8_t(l_v.B * c);
    }
    
    return Pixel(r, g , b);
}


struct Light {
    Light(const Vec3D& p, const float& i) : position(p), intensity(i) {}
    Vec3D position;
    float intensity;
 };

struct Ray 
{  
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
 
struct Cube 
{  
    Cube(const Vec3D &c_, const float& r_, const Material& m_) : center(c_), radius(r_), material(m_)  { 
        vertices.push_back(Vec3D(center.x - radius, center.y + radius, center.z - radius));
        vertices.push_back(Vec3D(center.x - radius, center.y + radius, center.z + radius));
        vertices.push_back(Vec3D(center.x + radius, center.y + radius, center.z + radius));
        vertices.push_back(Vec3D(center.x + radius, center.y + radius, center.z - radius));

        vertices.push_back(Vec3D(center.x - radius, center.y - radius, center.z - radius));
        vertices.push_back(Vec3D(center.x - radius, center.y - radius, center.z + radius));
        vertices.push_back(Vec3D(center.x + radius, center.y - radius, center.z + radius));
        vertices.push_back(Vec3D(center.x + radius, center.y - radius, center.z - radius));
        
        norms.push_back(cross(Vec3D(vertices[3] - vertices[0]), Vec3D(vertices[4] - vertices[0])).normalize()); //0423
        norms.push_back(cross(Vec3D(vertices[0] - vertices[1]), Vec3D(vertices[5] - vertices[1])).normalize()); //0145
        norms.push_back(cross(Vec3D(vertices[1] - vertices[2]), Vec3D(vertices[6] - vertices[2])).normalize()); //1265
        norms.push_back(cross(Vec3D(vertices[2] - vertices[3]), Vec3D(vertices[7] - vertices[3])).normalize()); //2376
        norms.push_back(cross(Vec3D(vertices[4] - vertices[5]), Vec3D(vertices[6] - vertices[5])).normalize()); //4567
        norms.push_back(cross(Vec3D(vertices[0] - vertices[3]), Vec3D(vertices[2] - vertices[3])).normalize()); //0123
        std::cout << norms[0] << std::endl;
        std::cout << norms[1] << std::endl;
        std::cout << norms[2] << std::endl;
        std::cout << norms[3] << std::endl;
        std::cout << norms[4] << std::endl;
        std::cout << norms[5] << std::endl;
        
        bounds[0] = vertices[4]; 
        bounds[1] = vertices[2]; 
    } 

    // Cube(const Vec3D &b0, const Vec3D &b1) { 
    //     bounds[0] = b0; 
    //     bounds[1] = b1; 
    //     center = Vec3D((b0.x + b1.x) / 2., (b0.y + b1.y) / 2., (b0.z + b1.z) / 2.);
    // } 
    std::vector<Vec3D> vertices;
    std::vector<Vec3D> norms;
    bool ray_intersect(const Ray&, float&) const;
    Vec3D normInPoint(const Vec3D&) const;
    Vec3D bounds[2]; 
    Vec3D center;
    float radius;
    Material material;
}; 

struct Sphere {
    Vec3D center;
    float radius;
    Material material;

    Sphere(const Vec3D &c_, const float &r_, const Material &m_) : center(c_), radius(r_), material(m_) {}
    bool ray_intersect(const Vec3D&, const Vec3D&, float&) const; 
};

bool Sphere::ray_intersect(const Vec3D &orig, const Vec3D &dir, float &t0) const 
{
        Vec3D L = center - orig;
        float tca = L*dir;
        float d2 = L*L - tca*tca;
        if (d2 > radius*radius) return false;
        float thc = sqrtf(radius*radius - d2);
        t0       = tca - thc;
        float t1 = tca + thc;
        if (t0 < 0) t0 = t1;
        if (t0 < 0) return false;
        return true;
    }

bool equals(const float& a, const float& b) {
    //std::cout << a << " VS " << b << std::endl;
    if (fabs(a - b) < 0.01) {
        //std::cout << a << " VS " << b << "==" << fabs(a - b) <<  " true" << std::endl;
        return true;
    }
    //std::cout << a << " VS " << b << "==" << fabs(a - b) <<  " false" << std::endl;
    return false;
}

Vec3D Cube::normInPoint(const Vec3D& point) const
{
    //int i1 = 0, i2 = 1, i3 = 2, i4 = 3;
    //float d, d1, d2, d3, d4;
    // d1 = std::sqrt(std::pow(vertices[0].x - point.x, 2) + std::pow(vertices[0].y - point.y, 2) + std::pow(vertices[0].z - point.z, 2));
    // d2 = std::sqrt(std::pow(vertices[1].x - point.x, 2) + std::pow(vertices[1].y - point.y, 2) + std::pow(vertices[1].z - point.z, 2));
    // d3 = std::sqrt(std::pow(vertices[2].x - point.x, 2) + std::pow(vertices[2].y - point.y, 2) + std::pow(vertices[2].z - point.z, 2));
    // d4 = std::sqrt(std::pow(vertices[3].x - point.x, 2) + std::pow(vertices[3].y - point.y, 2) + std::pow(vertices[3].z - point.z, 2));
    // for (int i = 4; i < 8; ++i) {
    //     d = std::sqrt(std::pow(vertices[i].x - point.x, 2) + std::pow(vertices[i].y - point.y, 2) + std::pow(vertices[i].z - point.z, 2));
    //     if (d < d1) {
    //         d4 = d3;
    //         d3 = d2;
    //         d2 = d1;
    //         d1 = d;
    //         i4 = i3;
    //         i3 = i2;
    //         i2 = i1;
    //         i1 = i;
    //     } else if (d < d2) {
    //         d4 = d3;
    //         d3 = d2;
    //         d2 = d;
    //         i4 = i3;
    //         i3 = i2;
    //         i2 = i;
    //     } else if (d < d3) {
    //         d4 = d3;
    //         d3 = d;
    //         i4 = i3;
    //         i3 = i;
    //     } else if (d < d4) {
    //         d4 = d;
    //         i4 = i;
    //     }
    // }
    
    //std::cout << std::endl << "This " << point << std::endl;
    // std::cout << vertices[0].z << std::endl;
    // std::cout << vertices[3].z << std::endl;
    // std::cout << vertices[4].z << std::endl;
    // std::cout << vertices[7].z << std::endl;

    if (equals(point.z, vertices[0].z) && equals(point.z, vertices[3].z) && 
            equals(point.z, vertices[4].z) && equals(point.z, vertices[7].z)) {
        return norms[0];
    }
    if (equals(point.x, vertices[0].x) && equals(point.x, vertices[1].x) && 
            equals(point.x, vertices[4].x) && equals(point.x, vertices[5].x)) {
        return norms[1];
    }
    //std::cout << "============" << std::endl;
    if (equals(point.z, vertices[1].z) && equals(point.z, vertices[2].z) && 
            equals(point.z, vertices[5].z) && equals(point.z, vertices[6].z)) {
        return norms[2];
    }
    if (equals(point.x, vertices[2].x) && equals(point.x, vertices[3].x) && 
            equals(point.x, vertices[7].x) && equals(point.x, vertices[6].x)) {
        return norms[3];
    }
    if (equals(point.y, vertices[4].y) && equals(point.y, vertices[5].y) && 
            equals(point.y, vertices[6].y) && equals(point.y, vertices[7].y)) {
        return norms[4];
    }
    if (equals(point.y, vertices[0].y) && equals(point.y, vertices[1].y) && 
            equals(point.y, vertices[2].y) && equals(point.y, vertices[3].y)) {
        return norms[5];
    }

    std::cout << "This " << point.z << std::endl;
    std::cout << vertices[0].z << std::endl;
    std::cout << vertices[3].z << std::endl;
    std::cout << vertices[4].z << std::endl;
    std::cout << vertices[7].z << std::endl;

    std::cout << norms[0] << std::endl;
    std::cout << norms[1] << std::endl;
    std::cout << norms[2] << std::endl;
    std::cout << norms[3] << std::endl;
    std::cout << norms[4] << std::endl;
    std::cout << norms[5] << std::endl;
    std::cout << norms[6] << std::endl;
}

bool Cube::ray_intersect(const Ray& r, float& t) const 
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
    if (fabs(x) < 0.0000000001) { x = 0.; } 
    float y = a.z * b.x - a.x * b.z;
    if (fabs(y) < 0.0000000001) { y = 0.; } 
    float z = a.x * b.y - a.y * b.x;
    if (fabs(z) < 0.0000000001) { z = 0.; } 
    return Vec3D(x, y, z);
}

float dot(const Vec3D& a, const Vec3D &b) 
{   
    return a.x * b.x + a.y * b.y + a.z + b.z;
}

#endif //__TEMPLATES_H__