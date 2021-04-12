#ifndef __TEMPLATES_H__
#define __TEMPLATES_H__

bool cube_intersect_bottom = false;

struct Pixel_final
{
    Pixel_final(){}
    Pixel_final(const char& r_, const char& g_, const char& b_) : R(r_), G(g_), B(b_) {}
    char R;
    char G;
    char B;
};

struct Vec4D
{
    Vec4D() {
        x = 0.;
        y = 0.;
        z = 0.;
        a = 0.;
    }
    Vec4D(const float& x_, const float& y_, const float& z_, const float& a_) : x(x_), y(y_), z(z_), a(a_) {}
    float x;
    float y;
    float z;
    float a;
};

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
    float length_squared() const;
    float x;
    float y;
    float z;
};

struct Pixel
{
    Pixel(){}
    Pixel(const float& r_, const float& g_, const float& b_) : R(r_), G(g_), B(b_) {}
    float R;
    float G;
    float B;
};

struct Material {
    Material(const float &r_, const Vec4D& a_, const Pixel& c_, const float &s_) : refractive_index(r_), albedo(a_), diffuse_color(c_), specular_exponent(s_) {}
    Material() : refractive_index(1), albedo(1, 0, 0, 0), diffuse_color(), specular_exponent() {}
    float refractive_index;
    Vec4D albedo;
    Pixel diffuse_color;
    float specular_exponent;
};

Vec3D cross (const Vec3D& a, const Vec3D& b); 
float dot(const Vec3D& a, const Vec3D &b);

Vec3D operator-(const Vec3D& l_v, const Vec3D& r_v)
{
    float x = l_v.x - r_v.x;
    float y = l_v.y - r_v.y;
    float z = l_v.z - r_v.z;
    return Vec3D(x, y, z);
}

Vec3D operator-(const Vec3D l_v)
{
    float x = -l_v.x;
    float y = -l_v.y;
    float z = -l_v.z;
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

float Vec3D::length_squared() const {
    return x * x + y * y + z * z;
}

float Vec3D::norm() {
    return std::sqrt(x * x + y * y + z * z);
}

Vec3D reflect(const Vec3D& I, const Vec3D& N) {
    return N * 2.f * dot(I, N) - I;
}

Vec3D Vec3D::normalize() {
    *this = (*this) * (1. / norm());
    return *this;
}

std::ostream& operator<<(std::ostream& out, const Vec3D& v) {
    out << v.x << " " << v.y << " " << v.z;
    return out ;
}

std::ostream& operator<<(std::ostream& out, const Pixel& v) {
    out << v.R << " " << v.G << " " << v.B;
    return out ;
}

Pixel operator*(const Pixel& l_v, float c) 
{
    float r, g, b;
    r = l_v.R * c;
    g = l_v.G * c;
    b = l_v.B * c;
    
    return Pixel(r, g, b);
}

Pixel operator+(const Pixel& l_v, const Pixel& r_v)
{
    float r, g, b;
    r = l_v.R + r_v.R;
    g = l_v.G + r_v.G;
    b = l_v.B + r_v.B;
    
    return Pixel(r, g, b);
}

Pixel operator+(const Pixel& l_v, float c)
{
    float r, g, b;
    r = l_v.R + c;
    g = l_v.G + c;
    b = l_v.B + c;
    
    return Pixel(r, g, b);
}

Pixel operator/(const Pixel& l_v, float c)
{
    float r, g, b;
    r = l_v.R / c;
    g = l_v.G / c;
    b = l_v.B / c;
    return Pixel(r, g, b);
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
    std::vector<Vec3D> vertices;
    std::vector<Vec3D> norms;
    bool ray_intersect(const Ray&, float&);
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
    bool ray_intersect(const Ray&, float&); 
    Vec3D normInPoint(const Vec3D&) const;
};

struct Surface 
{
    Surface(const float& y_, const float& size_, const Material &m1_, const Material &m2_) : y(y_), q_size(size_), material1(m1_), material2(m2_) {
        material = material1;
    }
    float y, q_size;
    Material material1, material2;
    Material material;
    bool ray_intersect(const Ray&, float&);
    Vec3D normInPoint(const Vec3D&) const;
};

bool Surface::ray_intersect(const Ray& ray, float &t)
{
    Vec3D N(0, -1, 0);
    Vec3D planePoint(0, y, 0);

    Vec3D diff = ray.orig - planePoint;
	float prod1 = dot(diff, N);
	float prod2 = dot(-ray.dir, N);
    //std::cout << ray.dir << "          " << N << std::endl;
    if (prod2 > 0.00001) {
        float prod3 = prod1 / prod2;
        t = prod3;
        Vec3D point = ray.orig + ray.dir * prod3;
        //std::cout << ray.orig + ray.dir * prod3 << std::endl;
        int a = 0;
        int n_x = int (point.x / q_size) + 1;
        int n_z = int (point.z / q_size) + 1;
        if (point.x < 0 && point.z >= 0 || point.x > 0 && point.z < 0) {
            a = 1;
        }
        if (((n_x + n_z + a) & 1) == 0) {
            material = material1;
        } else {
            material = material2;
        }
	    return true;
    }
    return false;
	
}



Vec3D Surface::normInPoint(const Vec3D& point) const
{
    return Vec3D(0, -1, 0);
}

Vec3D Sphere::normInPoint(const Vec3D& point) const
{
    return (point - center).normalize();
}

bool Sphere::ray_intersect(const Ray& ray, float &t) 
{
    Vec3D L = center - ray.orig;
    float tca = L*ray.dir;
    float d2 = dot(L, L) - tca*tca;
    if (d2 > radius*radius) return false;
    float thc = sqrtf(radius*radius - d2);
    t = tca - thc;
    float t1 = tca + thc;
    if (t < 0) t = t1;
    if (t < 0) return false;
    return true;
}

bool equals(const float& a, const float& b) {
    if (fabs(a - b) < 0.001) {
        return true;
    }
    return false;
}

Vec3D Cube::normInPoint(const Vec3D& point) const
{
    if (equals(point.z, vertices[0].z) && equals(point.z, vertices[3].z) && 
            equals(point.z, vertices[4].z) && equals(point.z, vertices[7].z)) {
        return norms[0];
    }
    if (equals(point.x, vertices[0].x) && equals(point.x, vertices[1].x) && 
            equals(point.x, vertices[4].x) && equals(point.x, vertices[5].x)) {
        return norms[1];
    }
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

    std::cout << "This " << point << std::endl;
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

bool Cube::ray_intersect(const Ray& r, float& t) 
{ 
    cube_intersect_bottom = false;
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
    //std::cout << tymin << " " << tymax << std::endl;
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
    //std::cout << r.orig + r.dir * t << std::endl;
    if (fabs((r.orig + r.dir * t).y - 4) < 0.01) {
        cube_intersect_bottom = true;
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
    return a.x * b.x + a.y * b.y + a.z * b.z;
}



struct Cone
{
    Cone(){}
    Cone(const Vec3D& c_, const float& r_, const float& h_, const Material& m1_, const Material& m2_) : material1(m1_), material2(m2_), center(c_), radius(r_), height(h_) {}
    bool ray_intersect(const Ray&, float&);
    Vec3D normInPoint(const Vec3D&) const;
    Vec3D center;
    float radius;
    float height;
    Material material;
    Material material1, material2;;
};

bool intersectPlane(const Ray& ray, const Vec3D &N, const Vec3D &p0, float &t) 
{ 
    Vec3D diff = ray.orig - p0;
	float prod1 = dot(diff, N);
	float prod2 = dot(-ray.dir, N);
    //std::cout << ray.dir << "          " << N << std::endl;
    if (prod2 > 0.00001) {
        float prod3 = prod1 / prod2;
        t = prod3;
        //Vec3D point = ray.orig + ray.dir * prod3;
        return true;
    } 
    return false; 
} 

bool intersectDisk(const Ray& ray, const Vec3D &p0, const float &radius, float& t) 
{ 
    Vec3D N(0, -1, 0);
    if (intersectPlane(ray, N, p0, t)) { 
        Vec3D p = ray.orig + ray.dir * t; 
        Vec3D v = p - p0; 
        float d2 = dot(v, v); 
        return (sqrtf(d2) <= radius); 
     } 
 
     return false; 
}

Vec3D Cone::normInPoint(const Vec3D& point) const
{
    // float r = sqrt((point.x-center.x)*(point.x-center.x) + (point.y-center.y)*(point.y-center.y));
    // return Vec3D(point.x-center.x, r*(radius / height), point.y-center.y).normalize();
    //return Vec3D(point.x - center.x, 0, point.z - center.z).normalize();
    //std::cout << point << std::endl;
    if (point.y < height - 0.001 + center.y && point.y > -height + 0.001 + center.y) {
        return Vec3D(point.x - center.x, 0, point.z - center.z).normalize();
    } else if (point.y == -height + center.y) {
        return Vec3D(0, -1, 0);
    }
    //return Vec3D(0, 1, 0);
}

bool Cone::ray_intersect(const Ray& ray, float& t)
{
    
    // float A = ray.orig.x - center.x;
    // float B = ray.orig.y - center.y;
    // float D = height - ray.orig.z + center.z;
    // float tan = (radius / height) * (radius / height);
    // float a = (ray.dir.x * ray.dir.x) + (ray.dir.y * ray.dir.y) - (tan*(ray.dir.z * ray.dir.z));
    // float b = (2*A*ray.dir.x) + (2*B*ray.dir.y) + (2*tan*D*ray.dir.z);
    // float c = (A*A) + (B*B) - (tan*(D*D));
    // float delta = b*b - 4*(a*c);
	// if(fabs(delta) < 0.001) return false; 
    // if(delta < 0.0) return false;
    // float t1 = (-b - sqrt(delta))/(2*a);
    // float t2 = (-b + sqrt(delta))/(2*a);
    // if (t1>t2) t = t2;
    // else t = t1;
    // float r = ray.orig.z + t*ray.dir.z;
    // if ((r > center.z) and (r < center.z + height)) return true;
    // else return false;

    // float a = (ray.dir.x * ray.dir.x) + (ray.dir.z * ray.dir.z);
    // float b = 2 * (ray.dir.x * (ray.orig.x - center.x) + ray.dir.z * (ray.orig.z - center.z));
    // float c = (ray.orig.x - center.x) * (ray.orig.x - center.x) + (ray.orig.z - center.z) * (ray.orig.z - center.z) - radius * radius;
    // float delta = b*b - 4*(a*c);
	// if(fabs(delta) < 0.0001) return false; 
    // if(delta < 0.0) return false;
    // float t1 = (-b - std::sqrt(delta))/(2*a);
    // float t2 = (-b + std::sqrt(delta))/(2*a);
    // if (t1 > t2) t = t2;
    // else t = t1;
    // std::cout << t << std::endl;
    // float r = ray.orig.y + t * ray.dir.y;
    // if ((r <= center.y) && (r >= center.y - height))return true;
    // else return false;
    // определяет предполагаемое положение пересечения (начало луча - центр цилиндра)
    
    // float EPS = 0.001;
    // Vec3D intersectionPoint = ray.orig - center;
    // bool isBelongToCylinderBase = false;
    // float ts1 = (height - ray.orig.y + center.y) / ray.dir.y;
    // Vec3D point = intersectionPoint + ray.dir * ts1;
    // if (point.x * point.x + point.z * point.z < radius * radius) {
    //     isBelongToCylinderBase = true;
    // }
    // float ts2 = (-height - ray.orig.y + center.y) / ray.dir.y;
    // point = intersectionPoint + ray.dir * ts2;
    // if (point.x * point.x + point.z * point.z < radius * radius) {
    //     isBelongToCylinderBase = true;
    // }
    material = material1;
    bool isBelongToCylinderBase = intersectDisk(ray, Vec3D(center.x, center.y - height, center.z), radius, t);
    if (isBelongToCylinderBase) {
        //std::cout << "!" << std::endl;
        material = material2;
        return true;
    }

    // определяет принадлежность пересечения цилиндру путем вычисления корня уравнения
    float a = ray.dir.x * ray.dir.x + ray.dir.z * ray.dir.z;
    float b = ray.orig.x * ray.dir.x - ray.dir.x * center.x + ray.orig.z * ray.dir.z - ray.dir.z * center.z;
    float c = ray.orig.x * ray.orig.x + center.x * center.x + ray.orig.z * ray.orig.z + center.z * center.z -
            2 * (ray.orig.x * center.x + ray.orig.z * center.z) - radius * radius;

    float discriminant = b * b - a * c;

    // если дискриминант меньше нуля, то пересечения с лучом нет и  t - бесконечность 
    if (discriminant < 0)
    {
        t = std::numeric_limits<float>::max();
        return false;
    }

    // вычисляет корни уравнения      
    float t1 = (-b - std::sqrtf(discriminant)) / a;
    float t2 = (-b + std::sqrtf(discriminant)) / a;

    // если корень меньше нуля, то пересечения с лучом нет 
    t = t1;
    if (t1 < 0)
        t = t2;
    //
    if (fabs(ray.orig.y + t * ray.dir.y - center.y) < height)  {
        return t > 0;
    }
    return false;
}

struct WorldObjects
{
    std::vector<Sphere> spheres;
    std::vector<Light> lights;
    std::vector<Cube> cubes;
    std::vector<Surface> surfaces;
    std::vector<Cone> cones;
    //std::vector<float> refr_indxs;
};

#endif //__TEMPLATES_H__