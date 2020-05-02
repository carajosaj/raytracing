#ifndef PETER_SHIRLEY_PROJECT_CODE_CONSTANT_MEDIUM_H
#define PETER_SHIRLEY_PROJECT_CODE_CONSTANT_MEDIUM_H

#include "hitable.h"
#include "texture.h"
#include "material.h"
#include"drand.h"

// 体，恒量介质
class constant_medium : public hitable {
public:
    constant_medium(hitable* b, double d, texture* a) : boundary(b), density(d) { phase_function = new isotropic(a); }
    virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const;
    virtual bool bounding_box(double t0, double t1, aabb& box) const {
        return boundary->bounding_box(t0, t1, box);
    }
    hitable* boundary;
    double density;
    // 材质为各项异性的材质
    material* phase_function;
};

bool constant_medium::hit(const ray& r, double t_min, double t_max, hit_record& rec) const {
    hit_record rec1, rec2;
    if (boundary->hit(r, -FLT_MAX, FLT_MAX, rec1)) {
        if (boundary->hit(r, rec1.t + 0.0001, FLT_MAX, rec2)) {
            if (rec1.t < t_min)
                rec1.t = t_min;
            if (rec2.t > t_max)
                rec2.t = t_max;
            if (rec1.t >= rec2.t)
                return false;
            if (rec1.t < 0)
                rec1.t = 0;
            double distance_inside_boundary = (rec2.t - rec1.t) * r.direction().length();
            double hit_distance = -(1 / density) * log(random_double());
            if (hit_distance < distance_inside_boundary) {
                rec.t = rec1.t + hit_distance / r.direction().length();
                rec.p = r.point_at_parameter(rec.t);
                rec.normal = vec3(1, 0, 0);  // arbitrary
                rec.mat_ptr = phase_function;
                return true;
            }
        }
    }
    return false;
}

#endif //PETER_SHIRLEY_PROJECT_CODE_CONSTANT_MEDIUM_H