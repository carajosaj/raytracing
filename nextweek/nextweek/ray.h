#ifndef RAYH //define rayh to avoid it be called for many times
#define RAYH
#include "vec3.h"

class ray
{
public:
    ray() {}
    ray(const vec3& a, const vec3& b,double ti=0.0) { A = a; B = b; _time = ti; }
    vec3 origin() const { return  A; }
    vec3 direction() const { return  B; }
    vec3 point_at_parameter(double t) const { return A + t * B; }

    double time() const { return  _time; }
    vec3 A;
    vec3 B;
    // 光线的时间戳
    double _time;

};

#endif 