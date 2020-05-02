#ifndef MATERIALH
#define MATERIALH

#include "ray.h"
#include "hitable.h"
#include"texture.h"
class material {
public:
	virtual bool scatter(
		const ray& r_in, const hit_record& rec, vec3& attenuation,
		ray& scattered) const = 0;
    virtual vec3 emitted(double u, double v, const vec3& p)const {
        return vec3(0, 0, 0);
    }
};
///////////////////////////////////////////////////////
class metal : public material
{
public:
    vec3 albedo;
    double fuzz;

    metal(const vec3& a, double f) : albedo(a)
    {
        fuzz = f < 1 ? f : 1;
    }

    virtual bool scatter(const ray& r_in, const hit_record& rec, vec3& attenuation, ray& scattered) const
    {
        vec3 v = unit_vector(r_in.direction());
        vec3 n = rec.normal;
        vec3 p = rec.p;
        vec3 r = reflect(v, n);
        vec3 offset = fuzz * random_in_unit_sphere();
        scattered = ray(p, r + offset);

        attenuation = albedo;
        return (dot(scattered.direction(), rec.normal) > 0);
    }
};
///////////////////////////////////////////////
class dielectric : public material {
public:
    dielectric(double ri) : ref_idx(ri) {}
    virtual bool scatter(
        const ray& r_in, const hit_record& rec, vec3& attenuation,
        ray& scattered) const
    {
        vec3 outward_normal;
        vec3 reflected = reflect(r_in.direction(), rec.normal);
        double ni_over_nt;
        attenuation = vec3(1.0, 1.0, 1.0);
        vec3 refracted;
        double reflect_prob;
        double cosine;
        if (dot(r_in.direction(), rec.normal) > 0) {
            outward_normal = -rec.normal;
            ni_over_nt = ref_idx;
            cosine = ref_idx * dot(r_in.direction(), rec.normal)
                / r_in.direction().length();
        }
        else {
            outward_normal = rec.normal;
            ni_over_nt = 1.0 / ref_idx;
            cosine = -dot(r_in.direction(), rec.normal)
                / r_in.direction().length();
        }
        if (refract(r_in.direction(), outward_normal, ni_over_nt, refracted)) {
            reflect_prob = schlick(cosine, ref_idx);
        }
        else {
            reflect_prob = 1.0;
        }
        if (random_double() < reflect_prob) {
            scattered = ray(rec.p, reflected);
        }
        else {
            scattered = ray(rec.p, refracted);
        }
        return true;
    }

    double ref_idx;
};
//////////////////////////////////////////////////
class lambertian : public material {//继承
public:
    lambertian(const vec3& a) : albedo(a) {}
    virtual bool scatter(const ray& r_in, const hit_record& rec,
        vec3& attenuation, ray& scattered) const
    {
        vec3 target = rec.p + rec.normal + random_in_unit_sphere();
        scattered = ray(rec.p, target - rec.p);
        attenuation = albedo;
        return true;
    }

    vec3 albedo;
};
//////////////////////////////////////////////
class lambertian1 : public material {
public:
    lambertian1(texture *a) : albedo(a) {}
    virtual bool scatter(const ray& r_in, const hit_record& rec, vec3& attenuation, ray& scattered) const {
        vec3 target = rec.p + rec.normal + random_in_unit_sphere();
        scattered = ray(rec.p, target - rec.p, r_in.time());
        attenuation = albedo->value(rec.u, rec.v, rec.p);
        return true;
    }

    texture* albedo;
};

class lambertian2 : public material {
public:
    lambertian2(texture* a) : albedo(a) {}
    virtual bool scatter(const ray& r_in, const hit_record& rec, vec3& attenuation, ray& scattered) const {
        vec3 target = rec.p + rec.normal + random_in_unit_sphere();
        scattered = ray(rec.p, target - rec.p);
        attenuation = albedo->value(0, 0, rec.p);
        return true;
    }

    texture* albedo;
};

//灯光材料
class diffuse_light :public material {
public:
    diffuse_light(texture* a) : emit(a) {}

    virtual bool scatter(const ray& r_in, const hit_record& rec, vec3& attenuation, ray& scattered) const {
        return false;
    }

    virtual vec3 emitted(double u, double v, const vec3& p) const {
        return emit->value(u, v, p);
    }

    texture* emit;
};

// 各向异性材质
class isotropic : public material {
public:
    isotropic(texture* a) : albedo(a) {}
    virtual bool scatter(const ray& r_in, const hit_record& rec, vec3& attenuation, ray& scattered) const {
        scattered = ray(rec.p, random_in_unit_sphere());
        attenuation = albedo->value(rec.u, rec.v, rec.p);
        return true;
    }
    texture* albedo;
};
#endif // !MATERIALH
