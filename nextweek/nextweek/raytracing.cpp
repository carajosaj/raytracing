// nextweek.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <iostream>
#include<fstream>
#include"camera.h"
#include"drand.h"
#include"material.h"
#define M_PI 3.14159265358979323846
#include"constant_mediun.h"


/*vec3 color(const ray& r, hitable* world, int depth) {
    hit_record rec;
    if (world->hit(r, 0.001, FLT_MAX, rec)) {
        ray scattered;
        vec3 attenuation;
        if (depth < 50 && rec.mat_ptr->scatter(r, rec, attenuation, scattered)) {
            return attenuation * color(scattered, world, depth + 1);
        }
        else {
            return vec3(0, 0, 0);
        }
    }
    else {
        vec3 unit_direction = unit_vector(r.direction());
        double t = 0.5 * (unit_direction.y() + 1.0);
        return (1.0 - t) * vec3(1.0, 1.0, 1.0) + t * vec3(0.5, 0.7, 1.0);
    }
}*/

vec3 color(const ray& r, hitable* world, int depth) {
    hit_record rec;
    if (world->hit(r, 0.001, FLT_MAX, rec)) {
        // 散射后的光线
        ray scattered;
        // 衰减
        vec3 attenuation;
        // 记录自发光的颜色
        vec3 emitted = rec.mat_ptr->emitted(rec.u, rec.v, rec.p);
        if (depth < 50 && rec.mat_ptr->scatter(r, rec, attenuation, scattered)) {
            // 递归 衰减
            return emitted + attenuation * color(scattered, world, depth + 1);
        }
        else {
            return emitted;
        }
    }
    else {
        //        vec3 unit_direction = unit_vector(r.direction());
        //        double t = 0.5 * (unit_direction.y() + 1.0);
        //        return (1.0 - t) * vec3(1.0, 1.0, 1.0) + t * vec3(0.5, 0.7, 1.0);
        return vec3(0, 0, 0);
    }
}

hitable* random_scene() {
    int n = 500;
    hitable** list = new hitable * [n + 1];
    list[0] = new sphere(vec3(0, -1000, 0), 1000, new lambertian(vec3(0.5, 0.5, 0.5)));
    int i = 1;
    for (int a = -11; a < 11; a++) {
        for (int b = -11; b < 11; b++) {
            double choose_mat = random_double();
            vec3 center(a + 0.9 * random_double(), 0.2, b + 0.9 * random_double());
            if ((center - vec3(4, 0.2, 0)).length() > 0.9) {
                if (choose_mat < 0.8) {  // diffuse
                    list[i++] = new sphere(center, 0.2,
                        new lambertian(vec3(random_double() * random_double(),
                            random_double() * random_double(),
                            random_double() * random_double())
                        )
                    );
                }
                else if (choose_mat < 0.95) { // metal
                    list[i++] = new sphere(center, 0.2,
                        new metal(vec3(0.5 * (1 + random_double()),
                            0.5 * (1 + random_double()),
                            0.5 * (1 + random_double())),
                            0.5 * random_double()));
                }
                else {  // glass
                    list[i++] = new sphere(center, 0.2, new dielectric(1.5));
                }
            }
        }
    }

    list[i++] = new sphere(vec3(0, 1, 0), 1.0, new dielectric(1.5));
    list[i++] = new sphere(vec3(-4, 1, 0), 1.0, new lambertian(vec3(0.4, 0.2, 0.1)));
    list[i++] = new sphere(vec3(4, 1, 0), 1.0, new metal(vec3(0.7, 0.6, 0.5), 0.0));

    return new hitable_list(list, i);
}
hitable* random_scene1() {
    int n = 500;
    hitable** list = new hitable * [n + 1];
    texture* checker = new checker_texture(new constant_texture(vec3(0.2, 0.3, 0.1)),
        new constant_texture(vec3(0.9, 0.9, 0.9)));
    list[0] = new sphere(vec3(0, -1000, 0), 1000, new lambertian2(checker));
    int i = 1;
    for (int a = -11; a < 11; a++) {
        for (int b = -11; b < 11; b++) {
            double choose_mat = random_double();
            vec3 center(a + 0.9 * random_double(), 0.2, b + 0.9 * random_double());
            if ((center - vec3(4, 0.2, 0)).length() > 0.9) {
                if (choose_mat < 0.8) {  // diffuse
                    // 运动模糊的小球
                    list[i++] = new moving_sphere(center, center + vec3(0, 0.5 * random_double(), 0), 0.0, 1.0, 0.2,
                        new lambertian(vec3(random_double() * random_double(), random_double() * random_double(),
                            random_double() * random_double())));
                }
                else if (choose_mat < 0.95) { // metal
                    list[i++] = new sphere(center, 0.2,
                        new metal(vec3(0.5 * (1 + random_double()), 0.5 * (1 + random_double()),
                            0.5 * (1 + random_double())), 0.5 * random_double()));
                }
                else {  // glass
                    list[i++] = new sphere(center, 0.2, new dielectric(1.5));
                }
            }
        }
    }

    list[i++] = new sphere(vec3(0, 1, 0), 1.0, new dielectric(2.5));
    list[i++] = new sphere(vec3(-4, 1, 0), 1.0, new lambertian(vec3(0.4, 0.2, 0.1)));
    list[i++] = new sphere(vec3(4, 1, 0), 1.0, new metal(vec3(1, 1, 1), 0.0));

    return new hitable_list(list, i);
}

hitable* two_spheres() {
    texture* checker = new checker_texture(new constant_texture(vec3(0.2, 0.3, 0.1)), 
        new constant_texture(vec3(0.9, 0.9, 0.9)));
    int n = 2;
    hitable** list = new hitable * [n + 1];
    list[0] = new sphere(vec3(0, -5, 0), 10, new lambertian1(checker));
    list[1] = new sphere(vec3(0, 5, 0), 10, new lambertian1(checker));

    return new hitable_list(list, 2);
}

hitable* two_perlin_spheres()
{
    texture* pertext = new noise_texture();
    hitable** list = new hitable * [2];
    list[0] = new sphere(vec3(0, -1000, 0), 1000, new lambertian2(pertext));
    list[1] = new sphere(vec3(0, 2, 0), 2, new lambertian2(pertext));
    return new hitable_list(list, 2);
}

hitable* simple_light()
{
    texture* pertext = new noise_texture(4);
    /*texture* checker = new checker_texture(new constant_texture(vec3(0.2, 0.3, 0.1)),
        new constant_texture(vec3(0.9, 0.9, 0.9)));*/
    hitable** list = new hitable * [4];
    list[0] = new sphere(vec3(0, -1000, 0), 1000, new lambertian2(pertext));
    list[1] = new sphere(vec3(0, 1, 0), 1, new lambertian2(pertext));
    list[2] = new sphere(vec3(0, 7, 0), 2, new diffuse_light(new constant_texture(vec3(4, 4, 4))));
    list[3] = new xy_rect(3, 5, 1, 3, -2, new diffuse_light(new constant_texture(vec3(4, 4, 4))));
    return new hitable_list(list, 4);
}

hitable* cornell_box() {
    hitable** list = new hitable * [8];
    int i = 0;
    material* red = new lambertian1(new constant_texture(vec3(0.65, 0.05, 0.05)));
    material* white = new lambertian1(new constant_texture(vec3(0.73, 0.73, 0.73)));
    material* green = new lambertian1(new constant_texture(vec3(0.12, 0.45, 0.15)));
    material* light = new diffuse_light(new constant_texture(vec3(15, 15, 15)));
    list[i++] = new flip_normals(new yz_rect(0, 555, 0, 555, 555, green));
    list[i++] = new yz_rect(0, 555, 0, 555, 0, red);
    list[i++] = new xz_rect(213, 343, 227, 332, 554, light);
    list[i++] = new flip_normals(new xz_rect(0, 555, 0, 555, 555, white));
    list[i++] = new xz_rect(0, 555, 0, 555, 0, white);
    list[i++] = new flip_normals(new xy_rect(0, 555, 0, 555, 555, white));
    //list[i++] = new box(vec3(130, 0, 65), vec3(295, 165, 230), white);
    //list[i++] = new box(vec3(265, 0, 295), vec3(430, 330, 460), white);
    list[i++] = new translate(new rotate_y(new box(vec3(0, 0, 0), vec3(165, 165, 165), white), -18), vec3(130, 0, 65));
    list[i++] = new translate(new rotate_y(new box(vec3(0, 0, 0), vec3(165, 330, 165), white), 15), vec3(265, 0, 295));
    return new hitable_list(list, i);
}

hitable* cornell_smoke()
{
    hitable** list = new hitable * [8];
    int i = 0;
    material* red = new lambertian2(new constant_texture(vec3(0.65, 0.05, 0.05)));
    material* white = new lambertian2(new constant_texture(vec3(0.73, 0.73, 0.73)));
    material* green = new lambertian2(new constant_texture(vec3(0.12, 0.45, 0.15)));
    material* light = new diffuse_light(new constant_texture(vec3(4, 4, 4)));
    list[i++] = new flip_normals(new yz_rect(0, 555, 0, 555, 555, green));
    list[i++] = new yz_rect(0, 555, 0, 555, 0, red);
    list[i++] = new xz_rect(113, 443, 127, 432, 554, light);
    list[i++] = new flip_normals(new xz_rect(0, 555, 0, 555, 555, white));
    list[i++] = new xz_rect(0, 555, 0, 555, 0, white);
    list[i++] = new flip_normals(new xy_rect(0, 555, 0, 555, 555, white));
    hitable* b1 = new translate(new rotate_y(new box(vec3(0, 0, 0), vec3(165, 165, 165), white), -18), vec3(130, 0, 65));
    hitable* b2 = new translate(new rotate_y(new box(vec3(0, 0, 0), vec3(165, 330, 165), white), 15), vec3(265, 0, 295));
    list[i++] = new constant_medium(b1,0.01,new constant_texture(vec3(1.0,1.0,1.0)));
    list[i++] = new constant_medium(b2, 0.01, new constant_texture(vec3(0.0, 0.0, 0.0)));
    return new hitable_list(list, i);
}
hitable* final() {
    int nb = 10;
    hitable** list = new hitable * [3000];
    material* white = new lambertian2(new constant_texture(vec3(0.73, 0.73, 0.73)));
    material* ground = new lambertian2(new constant_texture(vec3(0.48, 0.83, 0.53)));
    int b = 0;
    int l = 0;
    for (int i = 0; i < nb; i++) {
        for (int j = 0; j < nb; j++) {
            double w = 100;
            double x0 = i * w;
            double z0 = j * w;
            double y0 = 0;
            double x1 = x0 + w;
            double y1 = 100 * (random_double() + 0.01);
            double z1 = z0 + w;
            
            list[l++] = new box(vec3(x0, y0, z0), vec3(x1, y1, z1), ground);
        }
    }
    material* light = new diffuse_light(new constant_texture(vec3(7, 7, 7)));
    list[l++] = new xz_rect(123, 423, 147, 412, 554, light);
    vec3 center(400, 400, 200);
    list[l++] = new moving_sphere(center, center + vec3(30, 0, 0), 0, 1, 50,
        new lambertian2(new constant_texture(vec3(0.7, 0.3, 0.1))));
    list[l++] = new sphere(vec3(260, 150, 45), 50, new dielectric(1.5));
    list[l++] = new sphere(vec3(0, 150, 145), 50, new metal(vec3(0.8, 0.8, 0.9), 10.0));
    hitable* boundary = new sphere(vec3(360, 150, 145), 70, new dielectric(1.5));
    list[l++] = boundary;
    list[l++] = new constant_medium(boundary, 0.2, new constant_texture(vec3(0.2, 0.4, 0.9)));
    boundary = new sphere(vec3(0, 0, 0), 5000, new dielectric(1.5));
    list[l++] = new constant_medium(boundary, 0.0001, new constant_texture(vec3(1.0, 1.0, 1.0)));
    texture* pertext = new noise_texture(0.1);
    list[l++] = new sphere(vec3(220, 280, 300), 80, new lambertian2(pertext));
    int ns = 1000;
    for (int j = 0; j < ns; j++) {
        list[l++] = new sphere(vec3(165 * random_double() - 100, 165 * random_double() + 270, 165 * random_double() + 395), 10, white);
    }
   
    return new hitable_list(list, l);
}
int main()
{
    int nx = 200;
    int ny = 200;
    int ns = 10;
    ofstream fout("E:/ray1.ppm");
    fout << "P3\n" << nx << " " << ny << "\n255\n";
    double R = cos(M_PI / 4);
    vec3 lookfrom(278, 278, -800);
    vec3 lookat(278, 278, 0);
    double dist_to_focus = 10.0;
    double aperture = 0.1;
    double vfov = 40.0;
    camera cam(lookfrom, lookat, vec3(0, 1, 0), vfov, double(nx) / double(ny), aperture, dist_to_focus, 0.0, 1.0);
    /*vec3 lookfrom(13, 2, 3);
    vec3 lookat(0, 0, 0);
    double dist_to_focus = 10.0;
    double aperture = 0.0;

    camera cam(lookfrom, lookat, vec3(0, 1, 0), 20,
        double(nx) / double(ny), aperture, dist_to_focus, 0.0, 1.0);*/
    hitable* world = final();
    for (int j = ny - 1; j >= 0; j--) {
        for (int i = 0; i < nx; i++) {
            vec3 col(0, 0, 0);
            for (int s = 0; s < ns; s++) {
                double u = double(i + random_double()) / double(nx);
                double v = double(j + random_double()) / double(ny);
                ray r = cam.get_ray(u, v);
                vec3 p = r.point_at_parameter(2.0);
                col += color(r, world, 0);
            }
            col /= double(ns);
            col = vec3(sqrt(col[0]), sqrt(col[1]), sqrt(col[2]));
            int ir = int(255.99 * col[0]);
            int ig = int(255.99 * col[1]);
            int ib = int(255.99 * col[2]);
            ir = ir > 255 ? 255 : ir;
            ig = ig > 255 ? 255 : ig;
            ib = ib > 255 ? 255 : ib;
            fout << ir << " " << ig << " " << ib << "\n";
        }
    }
    fout.close();
}

