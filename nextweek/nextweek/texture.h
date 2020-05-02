#ifndef PETER_SHIRLEY_PROJECT_CODE_TEXTURE_H
#define PETER_SHIRLEY_PROJECT_CODE_TEXTURE_H

#include "vec3.h"
#include"perlin.h"

class texture {
public:
    virtual vec3 value(double u, double v, const vec3& p) const = 0;
};

class constant_texture : public texture {
public:
    constant_texture() {}

    constant_texture(vec3 c) : color(c) {}

    virtual vec3 value(double u, double v, const vec3& p) const {
        return color;
    }

    vec3 color;
};

//棋盘纹理
class checker_texture :public texture
{
public:
    checker_texture() {}
    checker_texture(texture* t0, texture* t1) :even(t0), odd(t1) {}
    virtual vec3 value(double u, double v, const vec3& p)const {
        double sines = sin(10 * p.x()) * sin(10 * p.y()) * sin(10 * p.z());
        if (sines < 0)
            return odd->value(u, v, p);
        else
            return  even->value(u, v, p);
    }
    // 棋盘纹理的间隔颜色
    texture* odd;
    texture* even;
};
//噪音纹理
class noise_texture :public texture {
public:
    noise_texture() {}
    noise_texture(double sc) :scale(sc) {}
    virtual vec3 value(double u, double v, const vec3& p)const
    {
        // 加缩放和扰动后
        return vec3(1, 1, 1) * 0.5 * (1 + sin(scale * p.z() + 10 * noise.turb(p)));
        //return vec3(1, 1, 1) *   noise.turb(scale*p);
    }
    perlin noise;
    double scale;
};

//照片纹理
class image_texture : public texture {
public:
    image_texture() {}
    image_texture(unsigned char* pixels, int A, int B) : data(pixels), nx(A), ny(B) {}
    virtual vec3 value(double u, double v, const vec3& p) const;
    unsigned char* data;
    int nx, ny;
};

vec3 image_texture::value(double u, double v, const vec3& p) const {
    int i = (1 - u) * nx;
    int j = (1 - v) * ny - 0.001;
    if (i < 0) i = 0;
    if (j < 0) j = 0;
    if (i > nx - 1) i = nx - 1;
    if (j > ny - 1) j = ny - 1;
    double r = int(data[3 * i + 3 * nx * j]) / 255.0;
    double g = int(data[3 * i + 3 * nx * j + 1]) / 255.0;
    double b = int(data[3 * i + 3 * nx * j + 2]) / 255.0;
    return vec3(r, g, b);
}

#endif //PETER_SHIRLEY_PROJECT_CODE_TEXTURE_H