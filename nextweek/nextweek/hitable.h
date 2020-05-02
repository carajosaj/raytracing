#ifndef HITABLEH
#define HITABLEH
#include "ray.h"
#include"aabb.h"
class material;

void get_sphere_uv(const vec3& p, double& u, double& v) {
    double phi = atan2(p.z(), p.x());
    double theta = asin(p.y());
    u = 1 - (phi + M_PI) / (2 * M_PI);
    v = (theta + M_PI / 2) / M_PI;
}

struct hit_record//记录参数
{
    double t;
    double u;
    double v;
    vec3 p;
    vec3 normal;
    material* mat_ptr;
};

class hitable
{
public:
    
	virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const = 0;
    virtual bool bounding_box(double t0, double t1, aabb& box)const = 0;

};

class hitable_list : public hitable {//多物体碰撞函数
public:
    hitable_list() {}
    hitable_list(hitable** l, int n) { list = l; list_size = n; }
    virtual bool hit(
        const ray& r, double tmin, double tmax, hit_record& rec) const;
    virtual bool bounding_box(double t0, double t1, aabb& box)const;
    hitable** list;
    int list_size;
};

bool hitable_list::hit(
    const ray& r, double t_min, double t_max, hit_record& rec) const {

    hit_record temp_rec;
    bool hit_anything = false;
    double closest_so_far = t_max;
    for (int i = 0; i < list_size; i++) {
        if (list[i]->hit(r, t_min, closest_so_far, temp_rec)) {
            hit_anything = true;
            closest_so_far = temp_rec.t;
            rec = temp_rec;
        }
    }
    return hit_anything;
}
bool hitable_list::bounding_box(double t0, double t1, aabb& box) const {
    if (list_size < 1) return false;
    aabb temp_box;
    bool first_true = list[0]->bounding_box(t0, t1, temp_box);
    if (!first_true)
        return false;
    else
        box = temp_box;
    for (int i = 1; i < list_size; i++) {
        if (list[0]->bounding_box(t0, t1, temp_box)) {
            box = surrounding_box(box, temp_box);
        }
        else
            return false;
    }
    return true;
}
////////////////////////////////////////////////////////
class sphere : public hitable {//球体碰撞函数
public:
    sphere() {}
    sphere(vec3 cen, double r, material* m)
        : center(cen), radius(r), mat_ptr(m) {};
    virtual bool hit(const ray& r, double tmin, double tmax, hit_record& rec) const;
    virtual bool bounding_box(double t0, double t1, aabb& box)const ;
    vec3 center;
    double radius;
    material* mat_ptr; /* NEW */
};

bool sphere::hit(const ray& r, double t_min, double t_max, hit_record& rec) const {
    vec3 oc = r.origin() - center;
    double a = dot(r.direction(), r.direction());
    double b = dot(oc, r.direction());
    double c = dot(oc, oc) - radius * radius;
    double discriminant = b * b - a * c;
    if (discriminant > 0) {
        double temp = (-b - sqrt(discriminant)) / a;
        if (temp < t_max && temp > t_min) {
            rec.t = temp;
            rec.p = r.point_at_parameter(rec.t);
            rec.normal = (rec.p - center) / radius;
            rec.mat_ptr = mat_ptr; /* NEW */
            return true;
        }
        temp = (-b + sqrt(discriminant)) / a;
        if (temp < t_max && temp > t_min) {
            rec.t = temp;
            rec.p = r.point_at_parameter(rec.t);
            rec.normal = (rec.p - center) / radius;
            rec.mat_ptr = mat_ptr; /* NEW */
            return true;
        }
    }
    return false;
}
bool sphere::bounding_box(double t0, double t1, aabb& box) const {
    box = aabb(center - vec3(radius, radius, radius), center + vec3(radius, radius, radius));
    return true;
}
/////////////////////////////////////////////////////////
//移动球体
class moving_sphere :public  hitable
{
public:
    moving_sphere() {}
    moving_sphere(vec3 cen0, vec3 cen1, double t0, double t1, double r, material* m) :
        center0(cen0), center1(cen1), time0(t0), time1(t1), radius(r), mat_ptr(m) {};

    virtual  bool hit(const ray& r, double tmin, double tmax, hit_record& rec) const;
    virtual bool bounding_box(double t0, double t1, aabb& box)const;
    vec3 center(double time) const;
    vec3 center0, center1;
    double time0, time1;
    double radius;
    material* mat_ptr;
};

// 当前时间点，球心的位置
vec3 moving_sphere::center(double time) const {
    return center0 + ((time - time0) / (time1 - time0)) * (center1 - center0);
}
//补充hit函数
bool moving_sphere::hit(const ray& r, double t_min, double t_max, hit_record& rec)const
{
    // 修改之前的center为一个时间相关的位置
    vec3 oc = r.origin() - center(r.time());
    double a = dot(r.direction(), r.direction());
    double b = dot(oc, r.direction());
    double c = dot(oc, oc) - radius * radius;
    double discriminant = b * b - a * c;
    if (discriminant > 0) {
        double temp = (-b - sqrt(discriminant)) / a;
        if (temp < t_max && temp > t_min) {
            rec.t = temp;
            rec.p = r.point_at_parameter(rec.t);
            rec.normal = (rec.p - center(r.time())) / radius;
            rec.mat_ptr = mat_ptr;
            return true;
        }
        temp = (-b + sqrt(discriminant)) / a;
        if (temp < t_max && temp > t_min) {
            rec.t = temp;
            rec.p = r.point_at_parameter(rec.t);
            rec.normal = (rec.p - center(r.time())) / radius;
            rec.mat_ptr = mat_ptr;
            return true;
        }
    }
    return false;
}
bool moving_sphere::bounding_box(double t0, double t1, aabb& box)const
{
    return true;
}
//////////////////////////////////////////////////////
class bvh_node :public hitable
{
public:
    bvh_node() {}
    bvh_node(hitable** l, int n, double time0, double time1);
    virtual bool hit(const ray& r, double tmin, double tmax, hit_record& rec)const;
    virtual bool bounding_box(double t0, double t1, aabb& box) const;

    hitable* left;
    hitable* right;
    aabb box;
};

bool bvh_node::hit(const ray& r, double tmin, double tmax, hit_record& rec) const {
    if (box.hit(r, tmin, tmax))
    {
        hit_record left_rec, right_rec;
        bool hit_left = left->hit(r, tmin, tmax, left_rec);
        bool hit_right = right->hit(r, tmin, tmax, right_rec);
        if (hit_left && hit_right)           // 击中重叠部分
        {
            if (left_rec.t < right_rec.t)
                rec = left_rec;             // 击中左子树
            else
                rec = right_rec;            // 击中右子树
            return true;
        }
        else if (hit_left)
        {
            rec = left_rec;
            return  true;
        }
        else if (hit_right)
        {
            rec = right_rec;
            return true;
        }
        else
            return false;
    }
    else
        return false;                       // 未击中任何物体
}
bool bvh_node::bounding_box(double t0, double t1, aabb& b) const {
    b = box;
    return true;
}

int box_x_compare(const void* a, const void* b)
{
    aabb box_left, box_right;
    hitable* ah = *(hitable**)a;
    hitable* bh = *(hitable**)b;
    if (!ah->bounding_box(0, 0, box_left) || !bh->bounding_box(0, 0, box_right))
        std::cerr << "No bounding box in bvh_node constructor\n";
    if (box_left.min().x() - box_right.min().x() < 0.0)
        return  -1;
    else
        return 1;
}


int box_y_compare(const void* a, const void* b)
{
    aabb box_left, box_right;
    hitable* ah = *(hitable**)a;
    hitable* bh = *(hitable**)b;
    if (!ah->bounding_box(0, 0, box_left) || !bh->bounding_box(0, 0, box_right))
        std::cerr << "no bounding box in bvh_node constructor\n";
    if (box_left.min().y() - box_right.min().y() < 0.0)
        return -1;
    else
        return 1;
}
int box_z_compare(const void* a, const void* b)
{
    aabb box_left, box_right;
    hitable* ah = *(hitable**)a;
    hitable* bh = *(hitable**)b;
    if (!ah->bounding_box(0, 0, box_left) || !bh->bounding_box(0, 0, box_right))
        std::cerr << "no bounding box in bvh_node constructor\n";
    if (box_left.min().z() - box_right.min().z() < 0.0)
        return -1;
    else
        return 1;
}

bvh_node::bvh_node(hitable** l, int n, double time0, double time1)
{
    int axis = int(3 * random_double());
    if (axis == 0)
        qsort(l, n, sizeof(hitable*), box_x_compare);
    else if (axis == 1)
        qsort(l, n, sizeof(hitable*), box_y_compare);
    else
        qsort(l, n, sizeof(hitable*), box_z_compare);
    if (n == 1) {
        left = right = l[0];
    }
    else if (n == 2) {
        left = l[0];
        right = l[1];
    }
    else {
        left = new bvh_node(l, n / 2, time0, time1);
        right = new bvh_node(l + n / 2, n - n / 2, time0, time1);
    }
    aabb box_left, box_right;
    if (!left->bounding_box(time0, time1, box_left) || !right->bounding_box(time0, time1, box_right))
        std::cerr << "no bounding box in bvh_node constructor\n";
    box = surrounding_box(box_left, box_right);
}
////////////////////////////////////////////////////////////
// xy平面的矩形
class xy_rect : public hitable {
public:
    xy_rect() {}
    xy_rect(double _x0, double _x1, double _y0, double _y1, double _k, material* mat) : x0(_x0), x1(_x1), y0(_y0), y1(_y1), k(_k), mp(mat) {};
    virtual bool hit(const ray& r, double t0, double t1, hit_record& rec) const;
    virtual bool bounding_box(double t0, double t1, aabb& box) const {
        box = aabb(vec3(x0, y0, k - 0.0001), vec3(x1, y1, k + 0.0001));
        return true;
    }
    material* mp;
    double x0, x1, y0, y1, k;
};
// 是否击中，形参传了hit_record的引用。
bool xy_rect::hit(const ray& r, double t0, double t1, hit_record& rec) const {
    double t = (k - r.origin().z()) / r.direction().z();
    if (t < t0 || t > t1)
        return false;
    double x = r.origin().x() + t * r.direction().x();
    double y = r.origin().y() + t * r.direction().y();
    if (x < x0 || x > x1 || y < y0 || y > y1)
        return false;
    rec.u = (x - x0) / (x1 - x0);
    rec.v = (y - y0) / (y1 - y0);
    rec.t = t;
    rec.mat_ptr = mp;
    rec.p = r.point_at_parameter(t);
    rec.normal = vec3(0, 0, 1);
    return true;
}
class xz_rect : public hitable {
public:
    xz_rect() {}
    xz_rect(double _x0, double _x1, double _z0, double _z1, double _k, material* mat) : x0(_x0), x1(_x1), z0(_z0), z1(_z1), k(_k), mp(mat) {};
    virtual bool hit(const ray& r, double t0, double t1, hit_record& rec) const;
    virtual bool bounding_box(double t0, double t1, aabb& box) const {
        box = aabb(vec3(x0, k - 0.0001, z0), vec3(x1, k + 0.0001, z1));
        return true;
    }
    material* mp;
    double x0, x1, z0, z1, k;
};

class yz_rect : public hitable {
public:
    yz_rect() {}
    yz_rect(double _y0, double _y1, double _z0, double _z1, double _k, material* mat) : y0(_y0), y1(_y1), z0(_z0), z1(_z1), k(_k), mp(mat) {};
    virtual bool hit(const ray& r, double t0, double t1, hit_record& rec) const;
    virtual bool bounding_box(double t0, double t1, aabb& box) const {
        box = aabb(vec3(k - 0.0001, y0, z0), vec3(k + 0.0001, y1, z1));
        return true;
    }
    material* mp;
    double y0, y1, z0, z1, k;
};

bool xz_rect::hit(const ray& r, double t0, double t1, hit_record& rec) const {
    double t = (k - r.origin().y()) / r.direction().y();
    if (t < t0 || t > t1)
        return false;
    double x = r.origin().x() + t * r.direction().x();
    double z = r.origin().z() + t * r.direction().z();
    if (x < x0 || x > x1 || z < z0 || z > z1)
        return false;
    rec.u = (x - x0) / (x1 - x0);
    rec.v = (z - z0) / (z1 - z0);
    rec.t = t;
    rec.mat_ptr = mp;
    rec.p = r.point_at_parameter(t);
    rec.normal = vec3(0, 1, 0);
    return true;
}

bool yz_rect::hit(const ray& r, double t0, double t1, hit_record& rec) const {
    double t = (k - r.origin().x()) / r.direction().x();
    if (t < t0 || t > t1)
        return false;
    double y = r.origin().y() + t * r.direction().y();
    double z = r.origin().z() + t * r.direction().z();
    if (y < y0 || y > y1 || z < z0 || z > z1)
        return false;
    rec.u = (y - y0) / (y1 - y0);
    rec.v = (z - z0) / (z1 - z0);
    rec.t = t;
    rec.mat_ptr = mp;
    rec.p = r.point_at_parameter(t);
    rec.normal = vec3(1, 0, 0);
    return true;
}
//反转法向量
class flip_normals : public hitable {
public:
    flip_normals(hitable* p) : ptr(p) {}
    virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const {
        if (ptr->hit(r, t_min, t_max, rec)) {
            rec.normal = -rec.normal;
            return true;
        }
        else
            return false;
    }
    virtual bool bounding_box(double t0, double t1, aabb& box) const {
        return ptr->bounding_box(t0, t1, box);
    }
    hitable* ptr;
};

//盒子
class box : public hitable {
public:
    box() {}
    box(const vec3& p0, const vec3& p1, material* ptr);
    virtual bool hit(const ray& r, double t0, double t1, hit_record& rec) const;
    virtual bool bounding_box(double t0, double t1, aabb& box) const {
        box = aabb(pmin, pmax);
        return true;
    }
    vec3 pmin, pmax;
    hitable* list_ptr;
};

box::box(const vec3& p0, const vec3& p1, material* ptr) {
    pmin = p0;
    pmax = p1;
    hitable** list = new hitable * [6];
    list[0] = new xy_rect(p0.x(), p1.x(), p0.y(), p1.y(), p1.z(), ptr);
    list[1] = new flip_normals(new xy_rect(p0.x(), p1.x(), p0.y(), p1.y(), p0.z(), ptr));
    list[2] = new xz_rect(p0.x(), p1.x(), p0.z(), p1.z(), p1.y(), ptr);
    list[3] = new flip_normals(new xz_rect(p0.x(), p1.x(), p0.z(), p1.z(), p0.y(), ptr));
    list[4] = new yz_rect(p0.y(), p1.y(), p0.z(), p1.z(), p1.x(), ptr);
    list[5] = new flip_normals(new yz_rect(p0.y(), p1.y(), p0.z(), p1.z(), p0.x(), ptr));
    list_ptr = new hitable_list(list, 6);
}

bool box::hit(const ray& r, double t0, double t1, hit_record& rec) const {
    return list_ptr->hit(r, t0, t1, rec);
}
///////////////////////////////////////////////////////////////
// 用于instance的移动
class translate : public hitable {
public:
    translate(hitable* p, const vec3& displacement) : ptr(p), offset(displacement) {}
    virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const;
    virtual bool bounding_box(double t0, double t1, aabb& box) const;
    hitable* ptr;
    vec3 offset;    // vec3的偏移
};

bool translate::hit(const ray& r, double t_min, double t_max, hit_record& rec) const {
    ray moved_r(r.origin() - offset, r.direction(), r.time());
    if (ptr->hit(moved_r, t_min, t_max, rec)) {
        rec.p += offset;
        return true;
    }
    else
        return false;
}

bool translate::bounding_box(double t0, double t1, aabb& box) const {
    if (ptr->bounding_box(t0, t1, box)) {
        box = aabb(box.min() + offset, box.max() + offset);
        return true;
    }
    else
        return false;
}
//旋转
class rotate_y : public hitable {
public:
    rotate_y(hitable* p, double angle);
    virtual bool hit(const ray& r, double t_min, double t_max, hit_record& rec) const;
    virtual bool bounding_box(double t0, double t1, aabb& box) const {
        box = bbox; return hasbox;
    }
    hitable* ptr;
    double sin_theta;
    double cos_theta;
    bool hasbox;
    aabb bbox;
};

rotate_y::rotate_y(hitable* p, double angle) : ptr(p) {
    double radians = (M_PI / 180.) * angle;
    sin_theta = sin(radians);
    cos_theta = cos(radians);
    hasbox = ptr->bounding_box(0, 1, bbox);
    vec3 min(FLT_MAX, FLT_MAX, FLT_MAX);
    vec3 max(-FLT_MAX, -FLT_MAX, -FLT_MAX);
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
                double x = i * bbox.max().x() + (1 - i) * bbox.min().x();
                double y = j * bbox.max().y() + (1 - j) * bbox.min().y();
                double z = k * bbox.max().z() + (1 - k) * bbox.min().z();
                double newx = cos_theta * x + sin_theta * z;
                double newz = -sin_theta * x + cos_theta * z;
                vec3 tester(newx, y, newz);
                // 旋转之后重新计算bounding box
                for (int c = 0; c < 3; c++)
                {
                    if (tester[c] > max[c])
                        max[c] = tester[c];
                    if (tester[c] < min[c])
                        min[c] = tester[c];
                }
            }
        }
    }
    bbox = aabb(min, max);
}

bool rotate_y::hit(const ray& r, double t_min, double t_max, hit_record& rec) const {
    vec3 origin = r.origin();
    vec3 direction = r.direction();
    //实现原点到方向之间的
    origin[0] = cos_theta * r.origin()[0] - sin_theta * r.origin()[2];
    origin[2] = sin_theta * r.origin()[0] + cos_theta * r.origin()[2];
    direction[0] = cos_theta * r.direction()[0] - sin_theta * r.direction()[2];
    direction[2] = sin_theta * r.direction()[0] + cos_theta * r.direction()[2];
    ray rotated_r(origin, direction, r.time());
    if (ptr->hit(rotated_r, t_min, t_max, rec)) {
        vec3 p = rec.p;
        vec3 normal = rec.normal;
        // normal 也做相应的旋转，因为是绕y轴，所以改p[0]和p[2]
        p[0] = cos_theta * rec.p[0] + sin_theta * rec.p[2];
        p[2] = -sin_theta * rec.p[0] + cos_theta * rec.p[2];
        normal[0] = cos_theta * rec.normal[0] + sin_theta * rec.normal[2];
        normal[2] = -sin_theta * rec.normal[0] + cos_theta * rec.normal[2];
        rec.p = p;
        rec.normal = normal;
        return true;
    }
    else
        return false;
}
//////////////////////////////////////////////////
// 体，恒量介质


#endif // !HITABLEH
