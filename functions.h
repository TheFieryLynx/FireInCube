#ifndef __FUNCTIONS_H__
#define __FUNCTIONS_H__

Vec3D cross (const Vec3D& a, const Vec3D& b) {
    Vec3D c;
    c.x = a.y * b.z - a.z * b.y;
    c.y = a.z * b.x - a.x * b.z;
    c.z = a.x * b.y - a.y * b.x;
    return c;
}

float dot(const Vec3D& a, const Vec3D &b) {
    return a.x * b.x + a.y * b.y + a.z + b.z;
}

//float triangle_intersection(const Vec3D& orig, const Vec3D& dir, const Vec3D& v0,
                                //const Vec3D& v1, const Vec3D& v2);

#endif //__FUNCTIONS_H__