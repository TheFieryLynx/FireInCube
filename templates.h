#ifndef __TEMPLATES_H__
#define __TEMPLATES_H__
#include <cmath>

template<size_t dimention, class Type>
struct vec 
{
    vec() {
        for (int i = 0; i < dimention; ++i) {
            data[i] = Type();
        }
    }
    Type& operator[](const size_t i) {
        return data[i];
    }
    const Type& operator[](const size_t i) const {
        return data[i];
    }
private:
    Type data[dimention];
};

typedef vec<3, float> Vec3D;
typedef vec<3, uint8_t> Pixel;
typedef vec<3, int> Rotation;

template<class Type>
struct vec<3, Type> 
{
    vec() : x(Type()), y(Type()), z(Type()) {}
    vec(Type x_, Type y_, Type z_) : x(x_), y(y_), z(z_) {}
    float norm() {
        return std::sqrt(x * x + y * y + z * z);
    };
    Type& operator[](const size_t i) {
        if (i == 0) {
            return x;
        } else if (i == 1) {
            return y;
        } else {
            return z;
        }
    }
    const Type& operator[](const size_t i) const {
        if (i == 0) {
            return x;
        } else if (i == 1) {
            return y;
        } else {
            return z;
        }
    }
    vec<3,Type> & normalize(Type l=1) { *this = (*this)*(l/norm()); return *this; }
    Type x;
    Type y;
    Type z;
};

template<size_t dimention,typename Type> 
Type operator*(const vec<dimention,Type>& lhs, const vec<dimention,Type>& rhs) {
    Type ret = Type();
    for (size_t i=dimention; i--; ret+=lhs[i]*rhs[i]);
    return ret;
}

template<size_t dimention,typename Type>
vec<dimention,Type> operator+(vec<dimention,Type> lhs, const vec<dimention,Type>& rhs) {
    for (size_t i=dimention; i--; lhs[i]+=rhs[i]);
    return lhs;
}

template<size_t dimention,typename Type>
vec<dimention,Type> operator-(vec<dimention,Type> lhs, const vec<dimention,Type>& rhs) {
    for (size_t i=dimention; i--; lhs[i]-=rhs[i]);
    return lhs;
}

template<size_t dimention,typename Type,typename U> 
vec<dimention,Type> operator*(const vec<dimention,Type> &lhs, const U& rhs) {
    vec<dimention,Type> ret;
    for (size_t i=dimention; i--; ret[i]=lhs[i]*rhs);
    return ret;
}

template <size_t dimention, typename Type> 
std::ostream& operator<<(std::ostream& out, const vec<dimention,Type>& v) {
    for(unsigned int i=0; i<dimention; i++) {
        out << v[i] << " " ;
    }
    return out ;
}

#endif //__TEMPLATES_H__