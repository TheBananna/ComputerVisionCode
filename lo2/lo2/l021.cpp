//Nicholas Bonanno PD: 3
#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <math.h>
#include <cmath>
#include <cstring>
#include <immintrin.h>
#include <smmintrin.h>
#include <iomanip>

using namespace std;
using namespace std::chrono;
struct vec2d//double point struct, but I call it a vec2d as that's what it technically is
{
    double x, y;
    vec2d()
    {
        x = 0.;
        y = 0.;
    }
    vec2d(double x, double y)
    {
        this->x = x;
        this->y = y;
    }
    vec2d operator*(int d)
    {
        return vec2d(this->x * d, this->y * d);
    }
    vec2d operator*(double d)
    {
        return vec2d(this->x * d, this->y * d);
    }
    vec2d operator -(vec2d vec)
    {
        return { x - vec.x, y - vec.y };
    }
    vec2d operator +(vec2d vec)
    {
        return { this->x + vec.x, this->y + vec.y };
    }
    vec2d operator/(double d)
    {
        return vec2d(this->x / d, this->y / d);
    }
    double mag()
    {
        return sqrt(x * x + y * y);
    }
    double mag2()
    {
        return (x * x + y * y);
    }
    double dist(vec2d vec)
    {
        return sqrt((vec.x - x) * (vec.x - x) + (vec.y - y) * (vec.y - y));
    }
    vec2d absolute()
    {
        return { abs(x), abs(y) };
    }
    vec2d midpoint(vec2d vec)
    {
        return (vec + *this) / 2;
    }
    double slope(vec2d vec)
    {
        return (vec.y - y) / (vec.x - x);
    }
    vec2d norm()
    {
        double length = mag();
        return {x / length, y / length};
    }
};
vec2d* get_points()
{
    vec2d* points = new vec2d[3];
    for (size_t i = 0; i < 3; i++)
    {
        points[i] = { (double)rand() / RAND_MAX, (double)rand() / RAND_MAX };
    }
    return points;
}
struct pixel//pixel, I'm using byte to reduce memory usage. This brings a measureable speedup, ~10% on my laptop
{
    uint8_t r, g, b;
};
const pixel BLACK = { 0, 0, 0 };
const pixel WHITE = { 255, 255, 255 };
inline void plot_pixel(int y, int x, pixel** pixels, pixel color = BLACK)
{
    //if (y >= 0 && x >= 0 && y < 800 && x < 800)
    if (((y | x) & 0x80000000) == 0 && y < 800 && x < 800)//bit trickery to reduce comparisons, I'm checking to see if either x or y have their negative bit set
        pixels[y][x] = color;
}
void draw_pixels(pixel** pixels)//writes pixels to a ppm file
{
    std::ofstream file("./pixels.ppm");
    if (!file.is_open())//in case something odd happens
    {
        cout << "FILE NOT OPENING";
        exit(1);
    }
    file.write("P3 800 800 255 \n", 16);
    
    for (size_t i = 0; i < 800; i++)
    {
        for (size_t j = 0; j < 800; j++)
        {
            file << std::to_string(pixels[i][j].r) << " " << std::to_string(pixels[i][j].g) << " " << std::to_string(pixels[i][j].b) << " ";
        }
        file << "\n";
    }
    file.flush();
    file.close();
}
struct vec2i//2 element int vector
{
    int x, y;
    vec2i()
    {
        x = 0;
        y = 0;
    }
    vec2i(int x, int y)
    {
        this->x = x;
        this->y = y;
    }
    vec2i(double x, double y)
    {
        this->x = round(x);
        this->y = round(y);
    }
    vec2i(vec2d vec)
    {
        this->x = round(vec.x);
        this->y = round(vec.y);
    }
};
void draw_horizontal_line(vec2i start, vec2i end, pixel** pixels)
{
    if (end.x < start.x)//so that it only goes left to right
    {
        auto temp = end;
        end = start;
        start = temp;
    }
    int dx = end.x - start.x;
    int dy = end.y - start.y;
    int yi = 1;
    if (dy < 0)
    {
        yi = -1;
        dy = -dy;
    }
    int D = 2 * dy - dx;
    int y = start.y;
    for (int x = start.x; x < end.x; x++)
    {
        plot_pixel(y, x, pixels);
        if (D > 0)
        {
            y = y + yi;
            D = D + 2 * (dy - dx);
        }
        else
        {
            D = D + 2 * dy;
        }
    }
    plot_pixel(end.y, end.x, pixels);
}
void draw_vertical_line(vec2i start, vec2i end, pixel** pixels)
{
    if (end.y < start.y)//so that it only goes down to up
    {
        auto temp = end;
        end = start;
        start = temp;
    }
    int dx = end.x - start.x;
    int dy = end.y - start.y;
    int xi = 1;
    if (dx < 0)
    {
        xi = -1;
        dx = -dx;
    }
    int D = 2 * dx - dy;
    int x = start.x;
    int y = start.y;
    for (; y < end.y; y++)
    {
        plot_pixel(y, x, pixels);
        if (D > 0)
        {
            x = x + xi;
            D = D + 2 * (dx - dy);
        }
        else
        {
            D = D + 2 * dx;
        }
    }
    plot_pixel(end.y, end.x, pixels);
}
void draw_line(vec2i start, vec2i end, pixel** pixels)//fills in the end point specifically
{
    if (abs(end.y - start.y) <= abs(end.x - start.x))
        draw_horizontal_line(start, end, pixels);
    else
        draw_vertical_line(start, end, pixels);
}
void draw_circle_brehnsham(vec2i center, int radius, pixel** pixels)
{
    plot_pixel(center.y, center.x, pixels);
    int r2 = (radius) * (radius)+radius;//to more accurately match an ideal circle made with y = sqrt(r^2 - x^2). ~equal to (x+.5)^2
    int y = radius;
    int y2 = y * y;
    int xmax = ceil(radius * 0.70710678118654752440084436210485) + 1;//ceil and +1 to ensure no undershoot
    int x2 = 0;
    for (size_t x = 0; x < xmax; x++)
    {
        if ((x2 + y2) > r2)//attempts to stay as close as possible to the ideal circle while not going outside of it
        {
            y -= 1;
            y2 = y * y;//recalculate y2
        }
        x2 = (x + 1) * (x + 1);
        plot_pixel(center.y + y, center.x + x, pixels);
        plot_pixel(center.y - y, center.x + x, pixels);
        plot_pixel(center.y + y, center.x - x, pixels);
        plot_pixel(center.y - y, center.x - x, pixels);

        plot_pixel(center.y + x, center.x + y, pixels);
        plot_pixel(center.y - x, center.x + y, pixels);
        plot_pixel(center.y + x, center.x - y, pixels);
        plot_pixel(center.y - x, center.x - y, pixels);

    }
}
inline double sign(vec2d p1, vec2d p2, vec2d p3)
{
    return (p1.x - p3.x) * (p2.y - p3.y) - (p2.x - p3.x) * (p1.y - p3.y);
}
bool in_triangle(vec2d p1, vec2d p2, vec2d p3, vec2d p4)
{
    double d1, d2, d3;
    bool neg, pos;

    d1 = sign(p4, p1, p2);
    d2 = sign(p4, p2, p3);
    d3 = sign(p4, p3, p1);

    neg = (d1 < 0) || (d2 < 0) || (d3 < 0);
    pos = (d1 > 0) || (d2 > 0) || (d3 > 0);

    return !(neg && pos);
}
inline uint8_t random_byte()
{
    return (uint8_t)(((double)rand() / (double)RAND_MAX) * 255);
}
vec2d part1(vec2d* points)
{
    vec2d point4;
    do
    {
        point4 = { (double)rand() / RAND_MAX, (double)rand() / RAND_MAX };
        //points = get_points();
        //vec2d midpoints[3] = { points[1].midpoint(points[2]), points[0].midpoint(points[2]), points[1].midpoint(points[0]) };
        //double sides[3];
        //sides[0] = ((points[1] - points[2]).absolute()).mag();//gets three side lengths, side[n] is opposite of points[n]
        //sides[1] = ((points[0] - points[2]).absolute()).mag();
        //sides[2] = ((points[0] - points[1]).absolute()).mag();

        //point4 = midpoints[0] + (midpoints[0] - points[0]).norm() * .1;
    } while ((in_triangle(points[0], points[1], points[2], point4) || in_triangle(points[1], points[2], point4, points[0]) || in_triangle(points[1], point4, points[0], points[2]) || in_triangle(point4, points[2], points[0], points[1])));
    //plot_pixel(point4.y * 800, point4.x * 800, pixels, { random_byte(), random_byte(), random_byte() });


    std::ofstream file("./points.txt");
    if (!file.is_open())//in case something odd happens
    {
        cout << "FILE NOT OPENING";
        exit(1);
    }
    file << setprecision(18);
    for (size_t i = 0; i < 3; i++)
    {
        file << "(" << points[i].x << "," << points[i].y << ") , ";
    }
    file << "(" << point4.x << "," << point4.y << ")";
    file.flush();
    file.close();
    return point4;
}
int main()
{
    //high_resolution_clock::time_point t1 = high_resolution_clock::now();
    srand(time(0));    
    part1(get_points());


    


    //pixel** pixels = new pixel * [800];//inits pixel grid to white
    //for (size_t i = 0; i < 800; i++)
    //{
    //    pixels[i] = new pixel[800];
    //    memset(pixels[i], 255, 800 * sizeof(pixel));
    //}
    //for (size_t i = 0; i < 3; i++)
    //{
    //    draw_circle_brehnsham({ points[i].x * 800, points[i].y * 800 }, 5, pixels);
    //}
    //draw_circle_brehnsham({ point4.x * 800, point4.y * 800 }, 10, pixels);

    //draw_line(points[1] * 800, points[2] * 800, pixels);
    //draw_line(points[2] * 800, points[0] * 800, pixels);
    //draw_line(points[0] * 800, points[1] * 800, pixels);

    //draw_pixels(pixels);
    //high_resolution_clock::time_point t2 = high_resolution_clock::now();
    //duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
    //std::cout << time_span.count();


}

