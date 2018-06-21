#ifndef __AF_VECTOR_TYPE_HELPERS_CUH__
#define __AF_VECTOR_TYPE_HELPERS_CUH__

#include <vector_types.h>
#include <vector_functions.h>
#include <math_functions.h>

#define __VECTOR_FUNCTIONS_DECL__ static __inline__ __host__ __device__


// Forward declarations for bool3
struct bool3;
__VECTOR_FUNCTIONS_DECL__ bool3 make_bool3(bool x, bool y, bool z);



__VECTOR_FUNCTIONS_DECL__ float3 operator+(const float3& a, const float3& b)
{
    return make_float3(a.x + b.x, a.y + b.y, a.z + b.z);
}

__VECTOR_FUNCTIONS_DECL__ float3 operator+(const float3& a, const float& s)
{
    return make_float3(a.x + s, a.y + s, a.z + s);
}

__VECTOR_FUNCTIONS_DECL__ uint3 operator+(const uint3& a, const uint3& b)
{
    return make_uint3(a.x + b.x, a.y + b.y, a.z + b.z);
}

__VECTOR_FUNCTIONS_DECL__ int3 operator+(const uint3& a, const int3& b)
{
    return make_int3(a.x + b.x, a.y + b.y, a.z + b.z);
}

__VECTOR_FUNCTIONS_DECL__ float3 operator-(const float3& a, const float3& b)
{
    return make_float3(a.x - b.x, a.y - b.y, a.z - b.z);
}

__VECTOR_FUNCTIONS_DECL__ void operator+=(float3& lhs, const float3& rhs)
{
    lhs.x += rhs.x; lhs.y += rhs.y; lhs.z += rhs.z;
}

__VECTOR_FUNCTIONS_DECL__ void operator-=(float3& lhs, const float3& rhs)
{
    lhs.x -= rhs.x; lhs.y -= rhs.y; lhs.z -= rhs.z;
}

__VECTOR_FUNCTIONS_DECL__ float3 operator*(const float3& a, const float& s)
{
    return make_float3(a.x * s, a.y * s, a.z * s);
}

__VECTOR_FUNCTIONS_DECL__ float3 operator/(const float3& a, const float& s)
{
    return make_float3(a.x / s, a.y / s, a.z / s);
}

__VECTOR_FUNCTIONS_DECL__ float3 operator*(const float& s, const float3& a)
{
    return make_float3(a.x * s, a.y * s, a.z * s);
}

__VECTOR_FUNCTIONS_DECL__ float3 operator/(const float& s, const float3& a)
{
    return make_float3(a.x / s, a.y / s, a.z / s);
}

__VECTOR_FUNCTIONS_DECL__ float3 operator/(const float3& a, const float3& b)
{
    return make_float3(a.x / b.x, a.y / b.y, a.z / b.z);
}

__VECTOR_FUNCTIONS_DECL__ int3 floor_float3(const float3& a)
{
    return make_int3(int(a.x), int(a.y), int(a.z));
}

__VECTOR_FUNCTIONS_DECL__ float mag(const float3& a)
{
    return sqrtf(a.x*a.x + a.y*a.y + a.z*a.z);
}

__VECTOR_FUNCTIONS_DECL__ float mag_sqrd(const float3& a)
{
    return a.x*a.x + a.y*a.y + a.z*a.z;
}

__VECTOR_FUNCTIONS_DECL__ float3 zero_float3()
{
    return make_float3(0,0,0);
}

__VECTOR_FUNCTIONS_DECL__ uint3 zero_uint3()
{
    return make_uint3(0,0,0);
}

__VECTOR_FUNCTIONS_DECL__ int3 operator%(const int3& a, const uint3& b)
{
    return make_int3(a.x % b.x, a.y % b.y, a.z % b.z);
}

//
// SOME CUSTOM VECTOR TYPES
//

struct bool3
{
    bool x;
    bool y;
    bool z;

    __inline__ __host__ __device__
    bool all(){ return x && y && z; }

    __inline__ __host__ __device__
    bool none(){ return (!x) && (!y) && (!z); }

    __inline__ __host__ __device__
    bool any(){ return x || y || z; }
};

__VECTOR_FUNCTIONS_DECL__ bool3 make_bool3(bool x, bool y, bool z)
{
    bool3 b; b.x = x; b.y = y; b.z = z;
    return b;
}

__VECTOR_FUNCTIONS_DECL__ bool3 operator>(const float3& a, const float3& b)
{
    return make_bool3(a.x > b.x, a.y > b.y, a.z > b.z);
}

__VECTOR_FUNCTIONS_DECL__ bool3 operator<(const float3& a, const float3& b)
{
    return make_bool3(a.x < b.x, a.y < b.y, a.z < b.z);
}

__VECTOR_FUNCTIONS_DECL__ bool3 operator>=(const float3& a, const float3& b)
{
    return make_bool3(a.x >= b.x, a.y >= b.y, a.z >= b.z);
}

__VECTOR_FUNCTIONS_DECL__ bool3 operator<=(const float3& a, const float3& b)
{
    return make_bool3(a.x <= b.x, a.y <= b.y, a.z <= b.z);
}

__VECTOR_FUNCTIONS_DECL__ bool3 operator&&(const bool3& a, const bool3& b)
{
    return make_bool3(a.x && b.x, a.y && b.y, a.z && b.z);
}

__VECTOR_FUNCTIONS_DECL__ bool3 operator||(const bool3& a, const bool3& b)
{
    return make_bool3(a.x || b.x, a.y || b.y, a.z || b.z);
}


__VECTOR_FUNCTIONS_DECL__ bool3 in_range(const int3& a, const uint3& l, const uint3& u)
{
    return make_bool3(
        a.x < u.x && a.x >= l.x, 
        a.y < u.y && a.y >= l.y, 
        a.z < u.z && a.z >= l.z);
}

#endif