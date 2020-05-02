#ifndef CAMERAH
#define CAMERAH

#include "drand.h"
#include "ray.h"
#define M_PI 3.14159265358979323846

vec3 random_in_unit_disk() {
    vec3 p;
    do {
        p = 2.0 * vec3(random_double(), random_double(), 0) - vec3(1, 1, 0);
    } while (dot(p, p) >= 1.0);
    return p;
}

class camera
{
    vec3 origin;
    vec3 u, v, w;
    vec3 horizontal;
    vec3 vertical;
    vec3 lower_left_corner;
    double len_radius;
    // 增加开始时间和结束时间
    double time0, time1;

public:
    // 构造函数增加t0，t1
    camera(vec3 lookfrom, vec3 lookat, vec3 vup, double vfov, double aspect, double aperture, double focus_dist,
        double t0, double t1)
    {
        time0 = t0;
        time1 = t1;
        len_radius = aperture / 2;
        double theta = vfov * M_PI / 180;
        double half_height = tan(theta / 2);
        double half_width = aspect * half_height;
        origin = lookfrom;

        w = unit_vector(lookfrom - lookat);
        u = unit_vector(cross(vup, w));
        v = cross(w, u);

        lower_left_corner = origin - half_width * focus_dist * u - half_height * focus_dist * v - focus_dist * w;
        horizontal = 2 * half_width * focus_dist * u;
        vertical = 2 * half_height * focus_dist * v;
    }

    ray get_ray(double s, double t)
    {
        vec3 rd = len_radius * random_in_unit_disk();
        vec3 offset = u * rd.x() + v * rd.y();
        // 随机时间戳的光线
        double time = time0 + random_double() * (time1 - time0);
        return ray(origin + offset, lower_left_corner + s * horizontal + t * vertical - origin - offset, time);
    }

};
#endif