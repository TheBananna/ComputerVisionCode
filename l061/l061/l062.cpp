//Nicholas Bonanno PD: 3
#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <math.h>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <vector>
#include <list>
#include <algorithm>
#include <tuple>
#include <random>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <deque>
#include <stack>
#include <regex>
#include <sstream>
#include <thread>
#include <mutex>
//#include <concurrent_unordered_map.h>
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
#include <Windows.h>
void sleep(int duration)
{
    Sleep(5000);
}
#else
#include <unistd.h>
void sleep(int duration)
{
    sleep(5000);
}
#endif

using namespace std;
using namespace std::chrono;

int IMG_HEIGHT;
int IMG_WIDTH;
int N;
int THRESHOLD2;
int MAGNITUDE_THRESHOLD;
int MAX_CENTER_OFFSET = 40;

template <class Key, class Value, class Hasher>
class concurrent_unordered_map
{
private:
    mutex mut;
public:
    unordered_map<Key, Value, Hasher> map;
    Value get (Key& k)
    {
        const lock_guard<mutex> lock(mut);
        return map[k];
    }
    void insert (Key k, Value val)
    {
        const lock_guard<mutex> lock(mut);
        map[k] = val;
    }
    bool contains(Key k)
    {
        const lock_guard<mutex> lock(mut);
        return map.find(k) != map.end();
    }
    Value operator [](Key k)
    {
        const lock_guard<mutex> lock(mut);
        return map.find(k) != map.end();
    }
    size_t size()
    {
        return map.size();
    }
};
class vec2d//double point struct, but I call it a vec2d as that's what it technically is
{
private:
    double x, y;
public:
    double get_x() { return x; };
    double get_y() { return y; };
    double set_x(double x) { this->x = x; return x; };
    double set_y(double y) { this->y = y; return y; };


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
    bool operator ==(vec2d v2) const
    {
        return (x == v2.get_x()) && (y == v2.get_y());
    }

    size_t operator()(const vec2d& point) const
    {
        return ((uint64_t)point.x ^ (uint64_t)point.y) + ((uint64_t)(point.x));
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
        return { x / length, y / length };
    }

    static vec2d rand_vec()
    {
        return { (double)rand() / RAND_MAX, (double)rand() / RAND_MAX };
    }
    double dist2(vec2d vec)
    {
        double dx = vec.x - x;
        double dy = vec.y - y;
        return dx * dx + dy * dy;
    }
    static int compare_x(const void* v1, const void* v2)
    {
        vec2d* vec1 = (vec2d*)v1;
        vec2d* vec2 = (vec2d*)v2;

        return (vec1->x > vec2->x) - (vec1->x < vec2->x);
    }
    static int compare_y(const void* v1, const void* v2)
    {
        vec2d* vec1 = (vec2d*)v1;
        vec2d* vec2 = (vec2d*)v2;

        return (vec1->y > vec2->y) - (vec1->y < vec2->y);
    }
    bool operator <(const vec2d& v2)
    {
        return x < v2.x;
    }
    bool operator >(const vec2d& v2)
    {
        return x > v2.x;
    }
    void print()
    {
        cout << x << " : " << y << endl;
    }
    double x_dist(vec2d v)
    {
        return abs(x - v.x);
    }
};
struct pixel//pixel, I'm using byte to reduce memory usage. This brings a measureable speedup, ~10% on my laptop
{
public:
    uint8_t r, g, b;
};
const pixel BLACK = { 0, 0, 0 };
const pixel WHITE = { 255, 255, 255 };
const pixel RED = { 255, 0, 0 };
const pixel GREEN = { 0, 255, 0 };
const pixel PURPLE = { 255, 0, 255 };
const pixel YELLOW = { 255, 255, 0 };
const pixel BLUE = { 0, 0, 255 };
inline void plot_pixel(int y, int x, pixel** pixels, pixel color = BLACK)
{
    //if (y >= 0 && x >= 0 && y < 400 && x < 400)
    if (((y | x) & 0x40000000) == 0 && y < IMG_HEIGHT && x < IMG_WIDTH)//bit trickery to reduce comparisons, I'm checking to see if either x or y have their negative bit set
        pixels[y][x] = color;
}
inline void plot_pixel(int y, int x, int** pixels)//plots white
{
    //if (y >= 0 && x >= 0 && y < 400 && x < 400)
    if (((y | x) & 0x40000000) == 0 && y < IMG_HEIGHT && x < IMG_WIDTH)//bit trickery to reduce comparisons, I'm checking to see if either x or y have their negative bit set
        pixels[y][x] = 255;
}
inline void add_pixel(int y, int x, int** pixels)
{
    //if (y >= 0 && x >= 0 && y < 400 && x < 400)
    if (((y | x) & 0x40000000) == 0 && y < IMG_HEIGHT && x < IMG_WIDTH)//bit trickery to reduce comparisons, I'm checking to see if either x or y have their negative bit set
        pixels[y][x] += 1;
}
void draw_pixels(pixel** pixels, string name, string max = "255")//writes pixels to a ppm file
{
    std::ofstream file("./" + name);
    if (!file.is_open())//in case something odd happens
    {
        cout << "FILE NOT OPENING";
        exit(1);
    }
    file.write(("P3 " + to_string(IMG_WIDTH) + " " + to_string(IMG_HEIGHT) + " " + max + " \n").c_str(), 7 + to_string(IMG_HEIGHT).size() + to_string(IMG_WIDTH).size() + max.size());

    for (size_t i = 0; i < IMG_HEIGHT; i++)
    {
        for (size_t j = 0; j < IMG_WIDTH; j++)
        {
            file << std::to_string(pixels[i][j].r) << " " << std::to_string(pixels[i][j].g) << " " << std::to_string(pixels[i][j].b) << " ";
        }
        file << "\n";
    }
}
void draw_pixels(int* pixels, string name, string max = "255")//writes pixels to a ppm file
{
    std::ofstream file("./" + name);
    if (!file.is_open())//in case something odd happens
    {
        cout << "FILE NOT OPENING";
        exit(1);
    }
    file.write(("P3 " + to_string(IMG_WIDTH) + " " + to_string(IMG_HEIGHT) + " " + max + " \n").c_str(), 7 + to_string(IMG_HEIGHT).size() + to_string(IMG_WIDTH).size() + max.size());

    for (size_t i = 0; i < IMG_HEIGHT * IMG_WIDTH; i++)
    {
        file << std::to_string(pixels[i]) << " " << std::to_string(pixels[i]) << " " << std::to_string(pixels[i]) << " ";
    }
    file << "\n";

}
void draw_pixels(int** pixels, string name, string max = "255")//writes pixels to a ppm file
{
    std::ofstream file("./" + name);
    if (!file.is_open())//in case something odd happens
    {
        cout << "FILE NOT OPENING";
        exit(1);
    }
    file.write(("P3 " + to_string(IMG_WIDTH) + " " + to_string(IMG_HEIGHT) + " " + max + " \n").c_str(), 7 + to_string(IMG_HEIGHT).size() + to_string(IMG_WIDTH).size() + max.size());

    for (size_t i = 0; i < IMG_HEIGHT; i++)
    {
        for (size_t j = 0; j < IMG_WIDTH; j++)
        {
            file << std::to_string(pixels[i][j]) << " " << std::to_string(pixels[i][j]) << " " << std::to_string(pixels[i][j]) << " ";
        }
        file << "\n";
    }
}
class vec2i//2 element int vector
{
private:
    int x, y;
public:
    int get_x() { return x; };
    int get_y() { return y; };
    int set_x(int x) { this->x = x; return x; };
    int set_y(int y) { this->y = y; return y; };

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
        this->x = round(vec.get_x());
        this->y = round(vec.get_y());
    }
    bool operator == (const vec2i& vec) const
    {
        return vec.x == x && vec.y == y;
    }
};
void draw_horizontal_line(vec2i start, vec2i end, pixel** pixels, pixel color = BLACK)
{
    if (end.get_x() < start.get_x())//so that it only goes left to right
    {
        auto temp = end;
        end = start;
        start = temp;
    }
    int dx = end.get_x() - start.get_x();
    int dy = end.get_y() - start.get_y();
    int yi = 1;
    if (dy < 0)
    {
        yi = -1;
        dy = -dy;
    }
    int D = 2 * dy - dx;
    int y = start.get_y();
    for (int x = start.get_x(); x < end.get_x(); x++)
    {
        plot_pixel(y, x, pixels, color);
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
    plot_pixel(end.get_y(), end.get_x(), pixels, color);
}
void draw_vertical_line(vec2i start, vec2i end, pixel** pixels, pixel color = BLACK)
{
    if (end.get_y() < start.get_y())//so that it only goes down to up
    {
        auto temp = end;
        end = start;
        start = temp;
    }
    int dx = end.get_x() - start.get_x();
    int dy = end.get_y() - start.get_y();
    int xi = 1;
    if (dx < 0)
    {
        xi = -1;
        dx = -dx;
    }
    int D = 2 * dx - dy;
    int x = start.get_x();
    int y = start.get_y();
    for (; y < end.get_y(); y++)
    {
        plot_pixel(y, x, pixels, color);
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
    plot_pixel(end.get_y(), end.get_x(), pixels, color);
}
void draw_line(vec2i start, vec2i end, pixel** pixels, pixel color = BLACK)//fills in the end point specifically
{
    if (abs(end.get_y() - start.get_y()) <= abs(end.get_x() - start.get_x()))
        draw_horizontal_line(start, end, pixels, color);
    else
        draw_vertical_line(start, end, pixels, color);
}
//#include <intrin.h>
void draw_circle_brehnsham(vec2i center, double radius, pixel** pixels, pixel color = BLACK)
{
    //unsigned long long start = __rdtsc();
    //plot_pixel(center.get_y(), center.get_x(), pixels);
    int r2 = round((radius) * (radius)+radius);//to more accurately match an ideal circle made with y = sqrt(r^2 - x^2). ~equal to (x+.5)^2
    int y = radius;
    int y2 = round(radius * radius);
    int xmax = ceil(radius * 0.70710678118654752440084436210485) + 1;//ceil and +1 to ensure no undershoot
    int x2 = 0;
    for (size_t x = 0; x < xmax; x++)
    {

        if ((x2 + y2) > r2)//attempts to stay as close as possible to the ideal circle while not going outside of it
        {
            y -= 1;
            y2 = y * y;//recalculate y2
        }
        double x1 = x + 1;
        x2 = x1 * x1;

        plot_pixel(center.get_y() + y, center.get_x() + x, pixels, color);
        plot_pixel(center.get_y() - y, center.get_x() + x, pixels, color);
        plot_pixel(center.get_y() + y, center.get_x() - x, pixels, color);
        plot_pixel(center.get_y() - y, center.get_x() - x, pixels, color);

        plot_pixel(center.get_y() + x, center.get_x() + y, pixels, color);
        plot_pixel(center.get_y() - x, center.get_x() + y, pixels, color);
        plot_pixel(center.get_y() + x, center.get_x() - y, pixels, color);
        plot_pixel(center.get_y() - x, center.get_x() - y, pixels, color);

    }
    //cout << __rdtsc() - start << endl;
}
void draw_circle_brehnsham(vec2i center, double radius, int** pixels)
{
    //unsigned long long start = __rdtsc();
    //plot_pixel(center.get_y(), center.get_x(), pixels);
    int r2 = round((radius) * (radius)+radius);//to more accurately match an ideal circle made with y = sqrt(r^2 - x^2). ~equal to (x+.5)^2
    int y = radius;
    int y2 = round(radius * radius);
    int xmax = ceil(radius * 0.70710678118654752440084436210485) + 1;//ceil and +1 to ensure no undershoot
    int x2 = 0;
    for (size_t x = 0; x < xmax; x++)
    {

        if ((x2 + y2) > r2)//attempts to stay as close as possible to the ideal circle while not going outside of it
        {
            y -= 1;
            y2 = y * y;//recalculate y2
        }
        double x1 = x + 1;
        x2 = x1 * x1;

        plot_pixel(center.get_y() + y, center.get_x() + x, pixels);
        plot_pixel(center.get_y() - y, center.get_x() + x, pixels);
        plot_pixel(center.get_y() + y, center.get_x() - x, pixels);
        plot_pixel(center.get_y() - y, center.get_x() - x, pixels);

        plot_pixel(center.get_y() + x, center.get_x() + y, pixels);
        plot_pixel(center.get_y() - x, center.get_x() + y, pixels);
        plot_pixel(center.get_y() + x, center.get_x() - y, pixels);
        plot_pixel(center.get_y() - x, center.get_x() - y, pixels);

    }
    //cout << __rdtsc() - start << endl;
}
void add_circle_brehnsham(vec2i center, double radius, int** pixels)
{
    //unsigned long long start = __rdtsc();
    //plot_pixel(center.get_y(), center.get_x(), pixels);
    int r2 = round((radius) * (radius)+radius);//to more accurately match an ideal circle made with y = sqrt(r^2 - x^2). ~equal to (x+.5)^2 if one considers the center of pixels to be the vertiecies of the grid
    int y = radius;
    int y2 = round(radius * radius);
    int xmax = ceil(radius * 0.70710678118654752440084436210485) + 1;//ceil and +1 to ensure no undershoot
    int x2 = 0;
    for (size_t x = 0; x < xmax; x++)
    {

        if ((x2 + y2) > r2)//attempts to stay as close as possible to the ideal circle while not going outside of it
        {
            y -= 1;
            y2 = y * y;//recalculate y2
        }
        double x1 = x + 1;
        x2 = x1 * x1;

        add_pixel(center.get_y() + y, center.get_x() + x, pixels);
        add_pixel(center.get_y() - y, center.get_x() + x, pixels);
        add_pixel(center.get_y() + y, center.get_x() - x, pixels);
        add_pixel(center.get_y() - y, center.get_x() - x, pixels);

        add_pixel(center.get_y() + x, center.get_x() + y, pixels);
        add_pixel(center.get_y() - x, center.get_x() + y, pixels);
        add_pixel(center.get_y() + x, center.get_x() - y, pixels);
        add_pixel(center.get_y() - x, center.get_x() - y, pixels);

    }
    //cout << __rdtsc() - start << endl;
}
void draw_incircle(vec2d p0, vec2d p1, vec2d p2, pixel** pixels)
{
    double sides[3];
    sides[0] = ((p1 - p2).absolute()).mag();//gets three side lengths, side[n] is opposite of points[n]
    sides[1] = ((p0 - p2).absolute()).mag();
    sides[2] = ((p0 - p1).absolute()).mag();
    double a = sides[0], b = sides[1], c = sides[2];//to reduce clutter
    double s = (a + b + c) / 2;
    vec2d center = { (a * p0.get_x() + b * p1.get_x() + c * p2.get_x()) / (a + b + c), (a * p0.get_y() + b * p1.get_y() + c * p2.get_y()) / (a + b + c) };//uses wikipedia algorithm based off of Heron's formula 
    double radius = round(sqrt((s - a) * (s - b) * (s - c) / s) * IMG_WIDTH);//if you want these to work you better be using a square I'm not mapping a rectangle to a square for you
    draw_circle_brehnsham(center * IMG_WIDTH, radius, pixels);
}
vec2d draw_circumcircle(vec2d p0, vec2d p1, vec2d p2, pixel** pixels)
{
    vec2d intercept;//algorithm I got off wikipedia
    double D = 2 * (p0.get_x() * (p1.get_y() - p2.get_y()) + p1.get_x() * (p2.get_y() - p0.get_y()) + p2.get_x() * (p0.get_y() - p1.get_y()));
    intercept.set_x((1 / D) * (p0.mag2() * (p1.get_y() - p2.get_y()) + p1.mag2() * (p2.get_y() - p0.get_y()) + p2.mag2() * (p0.get_y() - p1.get_y())));
    intercept.set_y((1 / D) * (p0.mag2() * (p2.get_x() - p1.get_x()) + p1.mag2() * (p0.get_x() - p2.get_x()) + p2.mag2() * (p1.get_x() - p0.get_x())));
    double radius = round(intercept.dist(p0) * IMG_WIDTH);
    draw_circle_brehnsham(intercept * IMG_WIDTH, radius, pixels);
    return intercept;
}
inline double sign(vec2d p1, vec2d p2, vec2d p3)
{
    return (p1.get_x() - p3.get_x()) * (p2.get_y() - p3.get_y()) - (p2.get_x() - p3.get_x()) * (p1.get_y() - p3.get_y());
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
void draw_line_to_edge(vec2d p1, vec2d p2, pixel** pixels)
{
    double slope1 = p1.slope(p2);//uses y = mx + b to allow easy calculation of start and end points, any pixels drawn outside of the 400x400 square will be discarded
    double b = p1.get_y() - slope1 * p1.get_x();
    draw_line(vec2d(0, b), vec2d(IMG_WIDTH, slope1 * IMG_HEIGHT + b), pixels);//could be more efficient but it doesn't really matter
}
pixel** pixels;

class line
{
private:
    vec2d p1;
    vec2d p2;
    double angle;
    double length;

public:
    vec2d get_p1() { return p1; }
    vec2d get_p2() { return p2; }
    double get_length() { return length; }
    double get_angle() { return angle; }
    void set_length(double length)
    {
        this->length = length;
        p2 = p1 + (vec2d(cos(angle), sin(angle)) * length);
    }
    line(vec2d p1, vec2d p2)
    {
        this->p1 = p1;
        this->p2 = p2;
        length = p1.dist(p2);
        angle = atan2(p2.get_y() - p1.get_y(), p2.get_x() - p1.get_x());
    }
    line()
    {
        p1 = {};
        p2 = {};
        angle = 0;
        length = 0;
    }
    void draw(pixel** pixels)
    {
        draw_line_to_edge(p1 * 400, p2 * 400, pixels);
        //draw_circle_brehnsham(p1 * 400, 5, pixels);
        //draw_circle_brehnsham(p2 * 400, 5, pixels);

    }
    vec2d intersection(line l)
    {
        //need to change to y=mx+b
        double slope1 = tan(angle);
        double b1 = p1.get_y() - slope1 * p1.get_x();

        double slope2 = tan(l.get_angle());
        double b2 = l.get_p1().get_y() - slope2 * l.p1.get_x();

        double intersect_x = (b2 - b1) / (slope1 - slope2);
        return vec2d(intersect_x, slope1 * intersect_x + b1);
    }
    double y_at(double x)
    {
        return (x - p1.get_x()) * tan(angle);
    }
    double left_right(vec2d v)//does a point fall left or right of the line?
    {
        return (p2.get_x() - p1.get_x()) * (v.get_y() - p1.get_y()) - (p2.get_y() - p1.get_y()) * (v.get_x() - p1.get_x());//I'll define right as negative
    }
    const double SQRT2 = 1.4142135623730950488016887242097;
    line get_perpendicular(vec2d p)
    {
        line perp;
        perp.angle = angle + M_PI_2 * (p.get_y() < y_at(p.get_x()) ? 1 : -1);
        perp.length = 10;
        perp.p1 = p;
        perp.p2 = p + (vec2d(cos(perp.angle), sin(perp.angle)) * 10);

        vec2d intersect = intersection(perp);
        perp.p1 = intersect;
        perp.p2 = intersect + (vec2d(cos(perp.angle), sin(perp.angle)) * intersect.dist(p1));
        perp.length = intersect.dist(perp.p2);
        draw_line(perp.p1 * IMG_WIDTH, perp.p2 * IMG_WIDTH, pixels);
        //draw_circle_brehnsham(p * 400, 7, pixels, BLUE);
        draw_line(p1 * IMG_WIDTH, p2 * IMG_WIDTH, pixels);
        return perp;
    }
    double distance_to_point(vec2d v)//returns squared value since I only use it for comparisons
    {
        return abs((p2.get_x() - p1.get_x()) * (p1.get_y() - v.get_y()) - (p1.get_x() - v.get_x()) * (p2.get_y() - p1.get_y())) / length;
    }
};

struct point_hasher
{
    size_t operator()(vec2i point) const
    {
        return (12119 + (point.get_x())) * 12119 + (point.get_y());//magic number is bilionth prime number
    }
};
bool contains(vector<vec2d> points, vec2d v)
{
    for (size_t i = 0; i < points.size(); i++)
    {
        if (points[i] == v)
            return true;
    }
    return false;
}
void find_hull(vector<vec2d> points, line l, vector<vec2d>* convex)
{
    int furthest = 0;
    double largest = 0;
    if (points.size() == 0)
        return;
    if (points.size() == 1)
    {
        convex->push_back(points[0]);
        return;
    }
    for (int i = 0; i < points.size(); i++)
    {
        double distance = l.distance_to_point(points[i]);
        if (distance > largest)
        {
            furthest = i;
            largest = distance;
        }
    }


    vec2d c = points[furthest];
    convex->push_back(points[furthest]);

    points.erase(points.begin() + furthest);
    vector<vec2d> left;
    vector<vec2d> right;
    line pc(l.get_p1(), c);
    line cq(c, l.get_p2());

    for (int i = 0; i < points.size(); i++)
    {
        if (!in_triangle(l.get_p1(), l.get_p2(), c, points[i]))
        {
            if (pc.left_right(points[i]) < 0) {
                left.push_back(points[i]);
            }
            else if (cq.left_right(points[i]) < 0) {
                right.push_back(points[i]);
            }
        }
    }
    find_hull(left, pc, convex);
    find_hull(right, cq, convex);
}
bool angle_pred(vec2d v1, vec2d v2) { return line(vec2d(0.5, 0.5), v1).get_angle() < line(vec2d(0.5, 0.5), v2).get_angle(); }
int THRESHOLD = 0;
class mat3x3
{
private:
    union {
        int mat[3][3];
        int row1[3];
        int row2[3];
        int row3[3];
    };
public:
    mat3x3(int mat[3][3])
    {
        this->mat[0][0] = mat[0][0];
        this->mat[0][1] = mat[0][1];
        this->mat[0][2] = mat[0][2];
        this->mat[1][0] = mat[1][0];
        this->mat[1][1] = mat[1][1];
        this->mat[1][2] = mat[1][2];
        this->mat[2][0] = mat[2][0];
        this->mat[2][1] = mat[2][1];
        this->mat[2][2] = mat[2][2];
    }
    void kernel(int* pixels, int** kernel)// evaluate a kernel on an image
    {
        for (int y = 1; y < IMG_HEIGHT - 1; y++)//ignore edges and sobel will be initialized black
        {
            for (int x = 1; x < IMG_WIDTH - 1; x++)
            {
                int sum = 0;
                for (int i = -1; i < 2; i++)//-1, 0, 1
                {
                    for (int j = -1; j < 2; j++)
                    {
                        sum += (mat[i + 1][j + 1] * (int)pixels[(y + i) * IMG_WIDTH + x + j]);
                    }
                }
                int val = sum;
                //cout << val << endl;
                //val = min(val > THRESHOLD ? val : 0, 255);
                kernel[y][x] = val;
            }
        }

    }
    void kernel(int** pixels, int** kernel)// evaluate a kernel on an image
    {
        for (int y = 1; y < IMG_HEIGHT - 1; y++)//ignore edges and sobel will be initialized black
        {
            for (int x = 1; x < IMG_WIDTH - 1; x++)
            {
                int sum = 0;
                for (int i = -1; i < 2; i++)//-1, 0, 1
                {
                    for (int j = -1; j < 2; j++)
                    {
                        sum += (mat[i + 1][j + 1] * (int)pixels[y + i][x + j]);
                    }
                }
                int val = sum;
                //cout << val << endl;
                //val = min(val > THRESHOLD ? val : 0, 255);
                kernel[y][x] = val;
            }
        }

    }

};
class mat5x5
{
private:
    union {
        int mat[5][5];
        int row1[5];
        int row2[5];
        int row3[5];
    };
public:
    mat5x5(int mat[5][5])
    {
        this->mat[0][0] = mat[0][0];
        this->mat[0][1] = mat[0][1];
        this->mat[0][2] = mat[0][2];
        this->mat[0][3] = mat[0][3];
        this->mat[0][4] = mat[0][4];
        this->mat[1][0] = mat[1][0];
        this->mat[1][1] = mat[1][1];
        this->mat[1][2] = mat[1][2];
        this->mat[1][3] = mat[1][3];
        this->mat[1][4] = mat[1][4];
        this->mat[2][0] = mat[2][0];
        this->mat[2][1] = mat[2][1];
        this->mat[2][2] = mat[2][2];
        this->mat[2][3] = mat[2][3];
        this->mat[2][4] = mat[2][4];
        this->mat[3][0] = mat[3][0];
        this->mat[3][1] = mat[3][1];
        this->mat[3][2] = mat[3][2];
        this->mat[3][3] = mat[3][3];
        this->mat[3][4] = mat[3][4];
        this->mat[4][0] = mat[4][0];
        this->mat[4][1] = mat[4][1];
        this->mat[4][2] = mat[4][2];
        this->mat[4][3] = mat[4][3];
        this->mat[4][4] = mat[4][4];
    }
    void kernel(int* pixels, int** kernel, int divisor)// evaluate a kernel on an image
    {
        for (int y = 2; y < IMG_HEIGHT - 2; y++)//ignore edges and sobel will be initialized black
        {
            for (int x = 2; x < IMG_WIDTH - 2; x++)
            {
                int sum = 0;
                for (int i = -2; i < 3; i++)//-1, 0, 1
                {
                    for (int j = -2; j < 3; j++)
                    {
                        sum += (mat[i + 2][j + 2] * (int)pixels[(y + i) * IMG_WIDTH + x + j]);
                    }
                }
                int val = sum;
                //cout << val << endl;
                //val = min(val > THRESHOLD ? val : 0, 255);
                kernel[y][x] = val / divisor;
            }
        }

    }
    void kernel(int** pixels, int** kernel, int divisor)// evaluate a kernel on an image
    {
        for (int y = 2; y < IMG_HEIGHT - 2; y++)//ignore edges and sobel will be initialized black
        {
            for (int x = 2; x < IMG_WIDTH - 2; x++)
            {
                int sum = 0;
                for (int i = -2; i < 3; i++)//-1, 0, 1
                {
                    for (int j = -2; j < 3; j++)
                    {
                        sum += (mat[i + 2][j + 2] * (int)pixels[y + i][x + j]);
                    }
                }
                int val = sum;
                //cout << val << endl;
                //val = min(val > THRESHOLD ? val : 0, 255);
                kernel[y][x] = val / divisor;
            }
        }

    }
};
int** get_image()//returns a 0 initialized image if IMG_WIDTH x IMG_HEIGHT
{
    int** img = new int* [IMG_HEIGHT];
    for (size_t i = 0; i < IMG_HEIGHT; i++)
    {
        img[i] = new int[IMG_WIDTH];
        memset(img[i], 0, IMG_WIDTH * sizeof(int));
    }
    return img;
}
unordered_set<vec2i, point_hasher> points;//to see if points have already been checked
long times = 0;
void flood_fill(vec2i pos, int** thresh, int depth)
{
    if (depth > 1000) return;//anti stack overflow protection in case it starts snaking
    if (points.find(pos) != points.end())
        return;
    points.insert(pos);

    times++;
    for (int y = -1; y < 2; y++)
    {
        for (int x = -1; x < 2; x++)
        {
            if ((((pos.get_y() + y > 0) && (pos.get_x() + x > 0))) && (y < IMG_HEIGHT) && (x < IMG_WIDTH) && (thresh[pos.get_y() + y][pos.get_x() + x] > 0))
                flood_fill(vec2i(pos.get_x() + x, pos.get_y() + y), thresh, depth + 1);
        }
    }
}
unordered_set<vec2i, point_hasher> points_groups;//to see if points have already been checked
void flood_fill_groups(vec2i pos, int** thresh, unordered_set<vec2i, point_hasher>& out, int depth)//gets groups of white pixels
{
    if (depth > 1000) return;//anti stack overflow protection in case it starts snaking
    if (thresh[pos.get_y()][pos.get_x()] <= 0 || points_groups.find(pos) != points_groups.end())//ignore this point if already checked or not valid
        return;
    points_groups.insert(pos);
    out.insert(pos);
    for (int y = -1; y < 2; y++)
    {
        for (int x = -1; x < 2; x++)
        {
            if ((((pos.get_y() + y >= 0) && (pos.get_x() + x >= 0))) && (pos.get_y() + y < IMG_HEIGHT) && (pos.get_x() + x < IMG_WIDTH) && (thresh[pos.get_y() + y][pos.get_x() + x] > 0))
                flood_fill_groups(vec2i(pos.get_x() + x, pos.get_y() + y), thresh, out, depth + 1);
        }
    }
}
int** copy_image(int** img)
{
    int** copy = get_image();
    for (int y = 0; y < IMG_HEIGHT; y++)
    {
        for (int x = 0; x < IMG_WIDTH; x++)
        {
            copy[y][x] = img[y][x];
        }
    }
    return copy;
}
//S is the kernel size for the second step and defines the maximum feature size detectable, should be even
//T is a threshold percent
bool REVERSE_BINARIZATION = false;
int** image_binarization(int** img, int S, double T)//based on Adaptive Thresholding Using the Integral Image byt Derek Bradley and Gerhard Roth
{
    int** binarizaiton = get_image();
    int** integral_image = get_image();

    for (int x = 0; x < IMG_WIDTH; x++)
    {
        int sum = 0;
        for (int y = 0; y < IMG_HEIGHT; y++)
        {
            sum += img[y][x];
            if (x == 0)
                integral_image[y][x] = sum;
            else
                integral_image[y][x] = integral_image[y][x - 1] + sum;
        }
    }
    for (int x = S / 2 + 1; x < IMG_WIDTH - S / 2; x++)
    {
        for (int y = S / 2 + 1; y < IMG_HEIGHT - S / 2; y++)
        {
            int x1 = x - S / 2;
            int x2 = x + S / 2;
            int y1 = y - S / 2;
            int y2 = y + S / 2;
            int count = (x2 - x1) * (y2 - y1);//this might possibly be a vector cross product not sure though
            int sum = integral_image[y2][x2] - integral_image[y1 - 1][x2] - integral_image[y2][x1 - 1] + integral_image[y1 - 1][x1 - 1];
            if ((img[y][x] * count) <= (sum * T))
                binarizaiton[y][x] = REVERSE_BINARIZATION ? 255 : 0;
            else
                binarizaiton[y][x] = REVERSE_BINARIZATION ? 0 : 255;
        }
    }
    return binarizaiton;
}
vec2i angle_map[9] = { {1, 0},{1, 1},{0, 1},{-1, 1},{-1, 0},{-1, -1},{0, -1},{1, -1},{1, 0} };//maps angle numbers to coordinate shifts
int** threshold(int** img, int thresh)
{
    int** thresholded = get_image();
    for (size_t y = 0; y < IMG_HEIGHT; y++)
    {
        for (size_t x = 0; x < IMG_WIDTH; x++)
        {
            if (img[y][x] >= thresh)
                thresholded[y][x] = img[y][x];
        }
    }
    return thresholded;
}
int** threshold_binary(int** img, int thresh)
{
    int** thresholded = get_image();
    for (size_t y = 0; y < IMG_HEIGHT; y++)
    {
        for (size_t x = 0; x < IMG_WIDTH; x++)
        {
            if (img[y][x] >= thresh)
                thresholded[y][x] = 1;
        }
    }
    return thresholded;
}
void part3(int argc, char** argv)
{
    //works if the photo comes from the convert utility mentioned in class, not guranteed otherwise
#pragma region FILE&ARG_READ
    const char* filename = "image.ppm";

    if (argc > 1)//file -
    {
        THRESHOLD = stoi(argv[2]);
        THRESHOLD2 = stoi(argv[4]);
        filename = argv[6];
    }
    ifstream ppm;
    ppm.open(filename);
    if (!ppm)
    {
        cout << "File doesn't exist!";
        exit(1);
    }
    string line;
    getline(ppm, line);
    getline(ppm, line);
    IMG_WIDTH = stoi(line.substr(0, line.find(' ')));
    IMG_HEIGHT = stoi(line.substr(line.find(' ') + 1));
    getline(ppm, line);
    int* pixels = new int[IMG_HEIGHT * IMG_WIDTH];//black and white
    int count = 0;
    int linecount = 0;
    while (getline(ppm, line))
    {
        //cout << line << ":" << endl << line.size() << endl;
        linecount++;
        //line = " " + line;//padding by one character
        int pos = 0;
        int r_pos = 0;
        int g_pos = 0;
        int b_pos = 0;
        while (pos < line.size())
        {
            r_pos = pos;
            g_pos = line.find(' ', pos + 1) + 1;
            b_pos = line.find(' ', g_pos) + 1;
            int third = line.find(' ', b_pos + 1) + 1;
            pos = third;//start of next r_pos


            int val = (((int)stoi(line.substr(r_pos, g_pos - r_pos - 1)) + (int)stoi(line.substr(g_pos, b_pos - g_pos - 1)) + (int)stoi(line.substr(b_pos, third - b_pos - 1))) / 3);
            pixels[count] = val;
            count++;
        }
    }
    cout << "here";
    draw_pixels(pixels, "imageg.ppm");
    cout << "grayscale done!" << endl;
#pragma endregion FILE&ARG_READ

    //I no longer map the sobel output onto [0, 255]
#pragma region SOBEL
    //const 
    //int* pixels_gaus = new int[IMG_HEIGHT * IMG_WIDTH];//black and white gausianed if requested
    //for (int y = 1; y < IMG_HEIGHT - 1; y++)
    //{
    //    for (int x = 1; x < IMG_WIDTH - 1; x++)
    //    {
    //        pixels_gaus[y * IMG_WIDTH + IMG_HEIGHT] = ()
    //    }
    //}


    int** sobelh = get_image();
    int** sobelv = get_image();
    int** sobel = get_image();

    int SOBELH[3][3] = { {-1, 0, 1}, {-2, 0, 2}, {-1, 0, 1} };
    int SOBELV[3][3] = { {1, 2, 1}, {0, 0, 0}, {-1, -2, -1} };
    mat3x3 mat_sobel_horiz = mat3x3(SOBELH);
    mat3x3 mat_sobel_vert = mat3x3(SOBELV);

    mat_sobel_horiz.kernel(pixels, sobelh);
    mat_sobel_vert.kernel(pixels, sobelv);

    int max = 0;

    for (int y = 0; y < IMG_HEIGHT; y++)
    {
        for (int x = 0; x < IMG_WIDTH; x++)
        {
            int val = round(sqrt(sobelh[y][x] * sobelh[y][x] + sobelv[y][x] * sobelv[y][x]));//.25 helps stop clipping, .125 is the theoretical scaling value based on the maximum output of a sobel operator but that rarely happens obviously so this also works well
            sobel[y][x] = val;//min(val > THRESHOLD ? val : 0, 255);
            if (val > max)
                max = val;
        }
    }

    //double scale = 255. / (double)max;
    //for (int y = 0; y < IMG_HEIGHT; y++)
    //{
    //    for (int x = 0; x < IMG_WIDTH; x++)
    //    {
    //        sobel[y][x] = min(255., sobel[y][x] * scale);//min(val > THRESHOLD ? val : 0, 255);
    //    }
    //}
    //draw_pixels(sobel, "imagem.ppm");//sobel but no threshold
    int** sobel_copy = copy_image(sobel);
    cout << "sobel done!" << endl;
#pragma endregion SOBEL

#pragma region DOUBLE_THRESHOLD
    int** thresh = get_image();
    for (size_t y = 0; y < IMG_HEIGHT; y++)
    {
        for (size_t x = 0; x < IMG_WIDTH; x++)
        {
            if (sobel[y][x] > THRESHOLD2)
                thresh[y][x] = 2;
            else if (sobel[y][x] > THRESHOLD)
                thresh[y][x] = 1;
            //inited to 0 so no need for a set to that
        }
    }


    for (int y = 0; y < IMG_HEIGHT; y++)
    {
        for (int x = 0; x < IMG_WIDTH; x++)
        {
            if (thresh[y][x] == 2)
            {
                flood_fill(vec2i(x, y), thresh, 0);
            }
        }
    }

    for (int y = 0; y < IMG_HEIGHT; y++)
    {
        for (int x = 0; x < IMG_WIDTH; x++)
        {
            if (points.find(vec2i(x, y)) == points.end())
            {
                sobel[y][x] = 0;
            }
            else
            {
                sobel[y][x] = 1;
            }
        }
    }
    draw_pixels(sobel, "image1.ppm", "1");//double threshold
    cout << "double threshold done!" << endl;
#pragma endregion DOUBLE_THRESHOLD

#pragma region NON_MAXIMUM_SUPPRESSION
    int** passed = get_image();
    for (int y = 1; y < IMG_HEIGHT - 1; y++)
    {
        for (int x = 1; x < IMG_WIDTH - 1; x++)
        {
            int angle = round(atan2((double)sobelv[y][x], (double)sobelh[y][x]) * 180. / M_PI / 45.);//divisions of 8 so that the rounding works in my favor here
            if (angle < 0)
                angle += 8;//no negative angles
            vec2i shift = angle_map[angle];

            //check to see if the supposed edge angle is lying or not
            if (sobel_copy[y][x] >= sobel_copy[y + shift.get_y()][x + shift.get_x()] && sobel_copy[y][x] >= sobel_copy[y - shift.get_y()][x - shift.get_x()])//if this is false fail the pixel
                passed[y][x] = 1;
            else
                passed[y][x] = 0;
        }
    }
    int** non_max_supp = copy_image(sobel);
    for (int y = 1; y < IMG_HEIGHT - 1; y++)
    {
        for (int x = 1; x < IMG_WIDTH - 1; x++)
        {
            if (!passed[y][x])
            {
                non_max_supp[y][x] = 0;
                sobel_copy[y][x] = 0;
            }
            else if (non_max_supp[y][x])
                non_max_supp[y][x] = 1;
        }
    }
    draw_pixels(passed, "image2.ppm", "1");
    cout << "non-maximum suppression done!" << endl;
    draw_pixels(non_max_supp, "imagef.ppm", "1");
    cout << "combination done!" << endl;
#pragma endregion NON_MAXIMUM_SUPPRESSION
}
int** canny(int** pixels, int THRESHOLD, int THRESHOLD2)
{
    //works if the photo comes from the convert utility mentioned in class, not guranteed otherwise

    //I no longer map the sobel output onto [0, 255]
#pragma region SOBEL
    //const 
    //int* pixels_gaus = new int[IMG_HEIGHT * IMG_WIDTH];//black and white gausianed if requested
    //for (int y = 1; y < IMG_HEIGHT - 1; y++)
    //{
    //    for (int x = 1; x < IMG_WIDTH - 1; x++)
    //    {
    //        pixels_gaus[y * IMG_WIDTH + IMG_HEIGHT] = ()
    //    }
    //}


    int** sobelh = get_image();
    int** sobelv = get_image();
    int** sobel = get_image();

    int SOBELH[3][3] = { {-1, 0, 1}, {-2, 0, 2}, {-1, 0, 1} };
    int SOBELV[3][3] = { {1, 2, 1}, {0, 0, 0}, {-1, -2, -1} };
    mat3x3 mat_sobel_horiz = mat3x3(SOBELH);
    mat3x3 mat_sobel_vert = mat3x3(SOBELV);

    mat_sobel_horiz.kernel(pixels, sobelh);
    mat_sobel_vert.kernel(pixels, sobelv);

    int max = 0;

    for (int y = 0; y < IMG_HEIGHT; y++)
    {
        for (int x = 0; x < IMG_WIDTH; x++)
        {
            int val = round(sqrt(sobelh[y][x] * sobelh[y][x] + sobelv[y][x] * sobelv[y][x]));//.25 helps stop clipping, .125 is the theoretical scaling value based on the maximum output of a sobel operator but that rarely happens obviously so this also works well
            sobel[y][x] = val;//min(val > THRESHOLD ? val : 0, 255);
            if (val > max)
                max = val;
        }
    }

    //double scale = 255. / (double)max;
    //for (int y = 0; y < IMG_HEIGHT; y++)
    //{
    //    for (int x = 0; x < IMG_WIDTH; x++)
    //    {
    //        sobel[y][x] = min(255., sobel[y][x] * scale);//min(val > THRESHOLD ? val : 0, 255);
    //    }
    //}
    //draw_pixels(sobel, "imagem.ppm");//sobel but no threshold
    int** sobel_copy = copy_image(sobel);
    //cout << "sobel done!" << endl;
#pragma endregion SOBEL

#pragma region DOUBLE_THRESHOLD
    int** thresh = sobel;
    for (size_t y = 0; y < IMG_HEIGHT; y++)
    {
        for (size_t x = 0; x < IMG_WIDTH; x++)
        {
            if (sobel[y][x] > THRESHOLD2)
                thresh[y][x] = 2;
            else if (sobel[y][x] > THRESHOLD)
                thresh[y][x] = 1;
            //inited to 0 so no need for a set to that
        }
    }


    for (int y = 0; y < IMG_HEIGHT; y++)
    {
        for (int x = 0; x < IMG_WIDTH; x++)
        {
            if (thresh[y][x] == 2)
            {
                flood_fill(vec2i(x, y), thresh, 0);
            }
        }
    }

    for (int y = 0; y < IMG_HEIGHT; y++)
    {
        for (int x = 0; x < IMG_WIDTH; x++)
        {
            if (points.find(vec2i(x, y)) == points.end())
            {
                sobel[y][x] = 0;
            }
            else
            {
                sobel[y][x] = 1;
            }
        }
    }
    //draw_pixels(sobel, "image1.ppm", "1");//double threshold
    //cout << "double threshold done!" << endl;
#pragma endregion DOUBLE_THRESHOLD

#pragma region NON_MAXIMUM_SUPPRESSION
    int** passed = get_image();
    for (int y = 1; y < IMG_HEIGHT - 1; y++)
    {
        for (int x = 1; x < IMG_WIDTH - 1; x++)
        {
            int angle = round(atan2((double)sobelv[y][x], (double)sobelh[y][x]) * 180. / M_PI / 45.);//divisions of 8 so that the rounding works in my favor here
            if (angle < 0)
                angle += 8;//no negative angles
            vec2i shift = angle_map[angle];

            //check to see if the supposed edge angle is lying or not
            if (sobel_copy[y][x] >= sobel_copy[y + shift.get_y()][x + shift.get_x()] && sobel_copy[y][x] >= sobel_copy[y - shift.get_y()][x - shift.get_x()])//if this is false fail the pixel
                passed[y][x] = 1;
            else
                passed[y][x] = 0;
        }
    }
    int** non_max_supp = copy_image(sobel);
    for (int y = 1; y < IMG_HEIGHT - 1; y++)
    {
        for (int x = 1; x < IMG_WIDTH - 1; x++)
        {
            if (!passed[y][x])
            {
                non_max_supp[y][x] = 0;
                sobel_copy[y][x] = 0;
            }
            else if (non_max_supp[y][x])
                non_max_supp[y][x] = 1;
        }
    }
    //draw_pixels(passed, "image2.ppm", "1");
    //cout << "non-maximum suppression done!" << endl;
    //draw_pixels(non_max_supp, "imagef.ppm", "1");
    //cout << "combination done!" << endl;
#pragma endregion NON_MAXIMUM_SUPPRESSION
    return non_max_supp;
    
}
int** summed_circular_hough_transform(int** img)
{
    int** hough_space = new int* [IMG_HEIGHT];//height x width x radius possibilities
    for (size_t z = 0; z < IMG_HEIGHT; z++)
    {
        hough_space[z] = new int[IMG_WIDTH];
        memset(hough_space[z], 0, 4 * IMG_WIDTH);
    }
    int iter = (int)floor((double)IMG_HEIGHT / (double)thread::hardware_concurrency());
    vector<thread*> threads;

    for (int thrd_id = 0; thrd_id < thread::hardware_concurrency(); thrd_id++)
    {
        thread* t = new thread([=]() {
            for (int y = iter * thrd_id; y < iter * (thrd_id + 1); y++)
            {
                if (y % (IMG_HEIGHT / 20) == 0)
                    cout << "y = " << to_string(y) << endl;
                for (int x = 0; x < IMG_WIDTH; x++)
                {
                    if (img[y][x] <= 0)
                        continue;
                    for (int R = 0; R < 120; R++)//125 is approximately the max radius of the coins
                    {
                        add_circle_brehnsham({ x, y }, R, hough_space);
                    }
                }
            }
            });
        threads.push_back(t);
    }
    for (thread* t : threads)
    {
        t->join();
        delete t;
    }

    int** hough_space_int = get_image();
    for (size_t y = 0; y < IMG_HEIGHT; y++)
    {
        for (size_t x = 0; x < IMG_WIDTH; x++)
        {
            hough_space_int[y][x] = hough_space[y][x];
        }
    }
    return hough_space_int;
}

int MAX_OFFSET = 50;
//returns point that is settled on starting at this point
vec2i gradient_ascent(int** sobelh, int** sobelv, int** img, int y, int x)
{
    vec2i past_pos = { -1, -1 };
    if (img[y][x] == 0)
        return past_pos;
    
    for (size_t i = 0; i < 300; i++)
    {
        if ((((y < MAX_OFFSET) || (x < MAX_OFFSET))) || (y >= IMG_HEIGHT - (MAX_OFFSET + 1)) || (x >= IMG_WIDTH - (MAX_OFFSET + 1)))
            return { -1, -1 };//gradient check will go off of the screen so don't count it

        vec2i largest = { x, y };
        for (int x_diff = -MAX_OFFSET; x_diff <= MAX_OFFSET; x_diff++)
        {
            for (int y_diff = -MAX_OFFSET; y_diff <= MAX_OFFSET; y_diff++)
            {
                if (img[largest.get_y()][largest.get_x()] < img[y + y_diff][x + x_diff])
                    largest = { x + x_diff, y + y_diff };
            }
        }

        if (largest == vec2i(x, y))
            return largest;//done
        past_pos = { x, y };
        y = largest.get_y();
        x = largest.get_x();

        //int angle = round(atan2((double)sobelv[y][x], (double)sobelh[y][x]) * 180. / M_PI / 45.);//divisions of 8 so that the rounding works in my favor here
        //if (angle < 0)
        //    angle += 8;//no negative angles
        //vec2i shift = angle_map[angle];
        //if (img[y][x] <= img[y + shift.get_y()][x + shift.get_x()])//if the gradient is correct
        //{
        //    past_pos = { x, y };
        //    y += shift.get_y();
        //    x += shift.get_x();
        //}
        //else
        //{
        //    return past_pos;//going to assume this is a maximum otherwise
        //}
    }
    return {-1, -1};//failed to find a maximum within the alloted travel time
}
int count_circle_edges(vec2i center, double radius, int** edges)
{
    //unsigned long long start = __rdtsc();
    //plot_pixel(center.get_y(), center.get_x(), pixels);
    int r2 = round((radius) * (radius)+radius);//to more accurately match an ideal circle made with y = sqrt(r^2 - x^2). ~equal to (x+.5)^2
    int y = radius;
    int y2 = round(radius * radius);
    int xmax = ceil(radius * 0.70710678118654752440084436210485) + 1;//ceil and +1 to ensure no undershoot
    int x2 = 0;
    int edge_count = 0;
    for (size_t x = 0; x < xmax; x++)
    {

        if ((x2 + y2) > r2)//attempts to stay as close as possible to the ideal circle while not going outside of it
        {
            y -= 1;
            y2 = y * y;//recalculate y2
        }
        double x1 = x + 1;
        x2 = x1 * x1;

        edge_count += edges[center.get_y() + y][center.get_x() + x];
        edge_count += edges[center.get_y() - y][center.get_x() + x];
        edge_count += edges[center.get_y() + y][center.get_x() - x];
        edge_count += edges[center.get_y() - y][center.get_x() - x];

        edge_count += edges[center.get_y() + x][center.get_x() + y];
        edge_count += edges[center.get_y() - x][center.get_x() + y];
        edge_count += edges[center.get_y() + x][center.get_x() - y];
        edge_count += edges[center.get_y() - x][center.get_x() - y];
    }
    return edge_count;
    //cout << __rdtsc() - start << endl;
}

tuple<vec2i, int> get_center(vec2i center_guess, int** edges)
{
    vec2i largest_group_pick = { center_guess.get_x(), center_guess.get_y() };
    int largest_group_pick_count = 0;
    int largest_group_pick_radius = 0;
    for (int radius = 70; radius < 120; radius++)
    {
        int x = center_guess.get_x(), y = center_guess.get_y();
        vec2i largest = { x, y };
        int largest_edge_count = count_circle_edges(largest, radius, edges);
        int largest_radius = 0;
        for (int x_diff = -MAX_OFFSET; x_diff <= MAX_OFFSET; x_diff++)
        {
            for (int y_diff = -MAX_OFFSET; y_diff <= MAX_OFFSET; y_diff++)
            {
                int edge_count = count_circle_edges({ x + x_diff, y + y_diff }, radius, edges);
                if (edge_count > largest_edge_count)
                {
                    largest = { x + x_diff, y + y_diff };
                    largest_edge_count = edge_count;
                    largest_radius = radius;
                }
            }
        }
        if (largest_edge_count > largest_group_pick_count)
        {
            largest_group_pick = largest;
            largest_group_pick_count = largest_edge_count;
            largest_group_pick_radius = largest_radius;
        }
    }
    return { largest_group_pick , largest_group_pick_radius };
}
void draw_bad_filled_circle(pixel** pixels, vec2i pos, int radius, pixel color)
{
    for (double i = 0; i <= radius; i += .5)
    {
        draw_circle_brehnsham(pos, i, pixels, color);
    }
}

void part2(int argc, char** argv)
{
    const char* filename = "image.ppm";
    //works if the photo comes from the convert utility mentioned in class, not guranteed otherwise
    //also grayscales the image
#pragma region FILE&ARG_READ


    if (argc > 1)//file -
    {
        THRESHOLD = stoi(argv[2]);
        THRESHOLD2 = stoi(argv[4]);
        filename = argv[6];
        cout << "It seems you've decided to use command line arguments. As a heads up the only one that actually affects anything here is the filename as I don't use Canny edge detection but the other two are still required for this to not crash." << endl << "sleeping for 7 seconds" << endl;
        sleep(7000);
    }
    ifstream ppm;
    ppm.open(filename);
    if (!ppm)
    {
        cout << "File doesn't exist!";
        exit(1);
    }
    string line;
    getline(ppm, line);
    getline(ppm, line);

    cout << ":" << line.substr(0, line.find(' ')) << ":" << endl;
    cout << ":" << line.substr(line.find(' ') + 1, line.size() - line.substr(0, line.find(' ')).size() - 1) << ":" << endl;;
    IMG_WIDTH = atoi(line.substr(0, line.find(' ')).c_str());
    IMG_HEIGHT = atoi(line.substr(line.find(' ') + 1, line.size() - line.substr(0, line.find(' ')).size() - 1).c_str());
    getline(ppm, line);

    int* pixels = new int[IMG_HEIGHT * IMG_WIDTH];//black and white
    pixel* pixels_color = new pixel[IMG_HEIGHT * IMG_WIDTH];
    int count = 0;
    int linecount = 0;
    while (getline(ppm, line))
    {
        //cout << linecount << endl;
        linecount++;
        //line = " " + line;//padding by one character
        int pos = 0;
        int r_pos = 0;
        int g_pos = 0;
        int b_pos = 0;
        while (pos < line.size())
        {
            r_pos = pos;
            g_pos = line.find(' ', pos + 1) + 1;
            b_pos = line.find(' ', g_pos) + 1;
            int third = line.find(' ', b_pos + 1) + 1;
            pos = third;//start of next r_pos

            int val = (((int)stoi(line.substr(r_pos, g_pos - r_pos - 1)) + (int)stoi(line.substr(g_pos, b_pos - g_pos - 1)) + (int)stoi(line.substr(b_pos, third - b_pos - 1))) / 3);
            pixels_color[count] = { (uint8_t)stoi(line.substr(r_pos, g_pos - r_pos - 1)), (uint8_t)stoi(line.substr(g_pos, b_pos - g_pos - 1)), (uint8_t)stoi(line.substr(b_pos, third - b_pos - 1)) };
            pixels[count] = val;
            count++;
        }
    }
    ppm.close();
    draw_pixels(pixels, "imageg.ppm");

    auto pixels_color_2d = new pixel*[IMG_HEIGHT];
    for (size_t y = 0; y < IMG_HEIGHT; y++)
    {
        pixels_color_2d[y] = new pixel[IMG_WIDTH];
        for (size_t x = 0; x < IMG_WIDTH; x++)
        {
            pixels_color_2d[y][x] = pixels_color[y * IMG_WIDTH + x];
        }
    }
    cout << "grayscale done!" << endl;
#pragma endregion FILE&ARG_READ

#pragma region HOUGH_TRANSFORM
    
    int gaussian_filter[5][5] = { {1, 4, 7, 4, 1}, { 4, 16, 26, 16, 4 }, { 7, 26, 41, 26, 7 }, { 4, 16, 26, 16, 4 }, { 1, 4, 7, 4, 1 } };
    mat5x5 gaussian = mat5x5(gaussian_filter);
    int** blurred = get_image();
    gaussian.kernel(pixels, blurred, 1);//no scaling back down to try to prevent zero clipping


    size_t blur_max = 0;
    for (int y = 0; y < IMG_HEIGHT; y++)//this doesn't work right, find out why
    {
        for (int x = 0; x < IMG_WIDTH; x++)
        {
            if (blurred[y][x] > blur_max)
                blur_max = blurred[y][x];
        }
    }
    draw_pixels(pixels, "blurred.ppm", to_string(blur_max));

    int unsharp_filter[5][5] = { {-1, -4, -6, -4, -1}, {-4, -16, -24, -16, -4}, {-6, -24, 476, -24, -6}, {-4, -16, -24, -16, -4}, {-1, -4, -6, -4, -1} };//{ {0, -1, 0}, {-1, 5, -1}, {0, -1, 5} };
    mat5x5 unsharp = mat5x5(unsharp_filter);
    int** sharpened = get_image();
    unsharp.kernel(blurred, sharpened, 256);//kernel was scaled by 256 from a floating point kernel to get this so I'll divide afterwards

    size_t avg_s = 0;
    size_t max_s = 0;
    for (int y = 0; y < IMG_HEIGHT; y++)//this doesn't work right, find out why
    {
        for (int x = 0; x < IMG_WIDTH; x++)
        {
            avg_s += sharpened[y][x];
            if (sharpened[y][x] > max_s)
                max_s = sharpened[y][x];
        }
    }
    cout << "sharpened avg " << avg_s / (IMG_HEIGHT * IMG_WIDTH) <<  ", maximum is " << max_s << endl;
    //draw_pixels(sharpened, "blurred.ppm");
    //int** thresholded = threshold(blurred, 1);
    //I no longer map the sobel output onto [0, 255]
    //int** binarization = image_binarization(sharpened, 8, .8);
    int** binarization = image_binarization(sharpened, 8, .8);//cahnge to .9?
    draw_pixels(binarization, "binar_orig.ppm", "1");
    int** binar_copy = /*canny(sharpened, THRESHOLD, THRESHOLD2);*/copy_image(binarization);
    int** binar_blurred = get_image();
    gaussian.kernel(binarization, binar_blurred, 1);//no scaling back down to try to prevent zero clipping
    binarization = binar_blurred;//memory leak I know and I know there's a lot more

    binar_blurred = threshold_binary(binarization, 1);
    binarization = binar_blurred;//image_binarization(binar_blurred, 8, .8);


    vector<unordered_set<vec2i, point_hasher>*> groups;
    for (int y = 0; y < IMG_HEIGHT; y++)//this doesn't work right, find out why
    {
        for (int x = 0; x < IMG_WIDTH; x++)
        {
            unordered_set<vec2i, point_hasher>* group = new unordered_set<vec2i, point_hasher>;
            flood_fill_groups(vec2i(x, y), binarization, *group, 0);
            //for (vec2i pos : *group)
            //{
            //    if (binarization[pos.get_y()][pos.get_x()] == 0);
            //    cout << "wrong";
            //}
            if (group->size() > 400)
                groups.push_back(group);
            else
            {
                for (vec2i pos : *group)
                    binarization[pos.get_y()][pos.get_x()] = 0;
                delete group;
            }
        }
    }
    draw_pixels(binarization, "imagef.ppm", "1");
    points_groups.clear();
    vector<unordered_set<vec2i, point_hasher>*> groups2;
    for (int y = 0; y < IMG_HEIGHT; y++)//this doesn't work right, find out why
    {
        for (int x = 0; x < IMG_WIDTH; x++)
        {
            unordered_set<vec2i, point_hasher>* group = new unordered_set<vec2i, point_hasher>;
            flood_fill_groups(vec2i(x, y), binar_copy, *group, 0);
            //for (vec2i pos : *group)
            //{
            //    if (binarization[pos.get_y()][pos.get_x()] == 0);
            //    cout << "wrong";
            //}
            if (group->size() > 200)
                groups2.push_back(group);
            else
            {
                for (vec2i pos : *group)
                    binar_copy[pos.get_y()][pos.get_x()] = 0;
                delete group;
            }
        }
    }
    
    draw_pixels(binar_copy, "edges.ppm", "1");
    //draw_pixels(binarization, "binary.ppm", "1");

    /*gaussian.kernel(binarization, binar_blurred, 273);
    binarization = copy_image(binar_blurred);
    gaussian.kernel(binarization, binar_blurred, 273);
    binarization = copy_image(binar_blurred);
    gaussian.kernel(binarization, binar_blurred, 273);
    binarization = copy_image(binar_blurred);
    gaussian.kernel(binarization, binar_blurred, 273);
    binarization = copy_image(binar_blurred);
    gaussian.kernel(binarization, binar_blurred, 273);
    binarization = copy_image(binar_blurred);*/

    //binarization = threshold_binary(binarization, 1);
    draw_pixels(binarization, "blurred.ppm");

    cout << "starting Hough Transform" << endl;
    int** hough = summed_circular_hough_transform(binarization);
    int max_found = 0;
    for (size_t i = 0; i < IMG_HEIGHT; i++)
    {
        for (size_t j = 0; j < IMG_WIDTH; j++)
        {
            max_found = max(max_found, hough[i][j]);
        }
    }
    //for (size_t i = 0; i < IMG_HEIGHT; i++)
    //{
    //    for (size_t j = 0; j < IMG_WIDTH; j++)
    //    {
    //        if (false)//(binarization[i][j] > 0)
    //            hough[i][j] = 255;
    //        else
    //            hough[i][j] = round(hough[i][j] * (255. / max));
    //    }
    //}
    //hough = threshold()
    //hough = image_binarization(hough, 4, .7);
    
    draw_pixels(hough, "imagev.ppm", to_string(max_found));

    


#pragma endregion HOUGH_TRANSOFRM

#pragma region GRADIENT_ASCENT
    cout << "Starting Gradient Ascent!" << endl;
    int** sobelh = get_image();
    int** sobelv = get_image();

    int SOBELH[3][3] = { {-1, 0, 1}, {-2, 0, 2}, {-1, 0, 1} };
    int SOBELV[3][3] = { {1, 2, 1}, {0, 0, 0}, {-1, -2, -1} };
    mat3x3 mat_sobel_horiz = mat3x3(SOBELH);
    mat3x3 mat_sobel_vert = mat3x3(SOBELV);

    mat_sobel_horiz.kernel(hough, sobelh);
    mat_sobel_vert.kernel(hough, sobelv);



    int iter = (int)floor((double)IMG_HEIGHT / (double)thread::hardware_concurrency());
    vector<thread*> threads;
    concurrent_unordered_map<vec2i, int, point_hasher> circles;
    auto circles_pointer = &circles;
    for (int thrd_id = 0; thrd_id < thread::hardware_concurrency(); thrd_id++)
    {
        thread* t = new thread([=]() {
            for (int y = iter * thrd_id; y < iter * (thrd_id + 1); y++)
            {
                if (y % (IMG_HEIGHT / 20) == 0)
                    cout << "y = " << to_string(y) << endl;
                for (int x = 0; x < IMG_WIDTH; x++)
                {
                    vec2i end = gradient_ascent(sobelh, sobelv, hough, y, x);
                    if (end == vec2i(-1, -1) /* || (end.get_x() < 250) || (end.get_y() < 87) || (end.get_x() > 3800) || (end.get_y() > 3000)*/)
                        continue;
                    else if (circles_pointer->contains(end))
                        circles_pointer->insert(end, circles_pointer->get(end) + 1);
                    else
                        circles_pointer->insert(end, 1);
                }
            }
            });
        threads.push_back(t);
    }
    for (thread* t : threads)
    {
        t->join();
        delete t;
    }
    threads.clear();

    int** circle_pixels = get_image();
    cout << circles.size() << " circles to draw" << endl;
    int drawn = 1;
    double avg = 0;
    vector<tuple<vec2i, int>> radii;
    for (std::pair<vec2i, int> pair : circles.map)
    {
        if (pair.first.get_x() < 210 || pair.first.get_x() > (IMG_WIDTH - 210) || pair.first.get_y() < 210 || pair.first.get_y() > (IMG_HEIGHT - 210))
            continue;//this is needed as that's the max radius checked by my get_center method + the max distance from the center guess
        auto circle = get_center(pair.first, binar_copy);
        if (get<1>(circle) < 40)
            continue;//circle is too small to realistically be a coin
        avg += get<1>(pair);
        radii.push_back(circle);
    }
    avg /= radii.size();
    double diff_sum = 0;
    for (std::tuple<vec2i, int> pair : radii)
        diff_sum += pow(get<1>(pair) - avg, 2);
    diff_sum /= radii.size();
    double std_dev = sqrt(diff_sum);
    sort(radii.begin(), radii.end(), [](tuple<vec2i, int> p1, tuple<vec2i, int> p2) {
        return get<1>(p1) > get<1>(p2);
    });
    //now I need to figure out which coin is which now, I'm going to attempt some basic grouping
    pixel colors[] = {YELLOW, GREEN, PURPLE, RED, BLUE};//this will only work properly for the easy or medium image
    int ind = 0;
    int size_prev = get<1>(radii[0]);
    int coins[5] = { 0,0,0,0,0 };//penny, nickel, quarter
    for (std::tuple<vec2i, int> circle : radii)
    {
        if (size_prev - get<1>(circle) > 1)
            ind = min(ind+1, 4);
        size_prev = get<1>(circle);
        coins[ind]++;
        //draw_circle_brehnsham(pair.first, pair.second / 10, circle_pixels);
        //[pair.first.get_y()][pair.first.get_x()] = pair.second;
        for (double rad = get<1>(circle) - 5.; rad < get<1>(circle); rad+=.5)
        {
            draw_circle_brehnsham(get<0>(circle), rad, pixels_color_2d, colors[ind]);
        }
        cout << "circles drawn " << drawn++ << "!" << endl;
    }
    for (int y_diff = 0; y_diff < 5; y_diff++)
    {
        draw_line({ 210, 210 + y_diff }, { IMG_WIDTH - 210, 210 + y_diff }, pixels_color_2d);
        draw_line({ 210, IMG_HEIGHT - 210 + y_diff }, { IMG_WIDTH - 210, IMG_HEIGHT - 210 + y_diff }, pixels_color_2d);
    }
    for (int x_diff = 0; x_diff < 5; x_diff++)
    {
        draw_line({ 210 + x_diff, 210 }, { 210 + x_diff, IMG_HEIGHT - 210 }, pixels_color_2d);
        draw_line({ IMG_WIDTH - 210 + x_diff, 210 }, { IMG_WIDTH - 210 + x_diff, IMG_HEIGHT - 210 }, pixels_color_2d);
    }
    draw_pixels(pixels_color_2d, "coins.ppm");

    //draw_pixels(circle_pixels, "circles.ppm");
#pragma endregion GRADIENT_ASCENT

    ofstream output("results.txt");
    cout << "Images lacking Dimes and Dollar coins won't detect properly as the coin detection algorithm goes down the sizes" << endl;
    cout << "As a consequence, this will only detect coins properly in the hard image but at the same time will fail to detect the coins properly in that image as well" << endl;
    string coin_names[5] = { "Dollar", "Quarter", "Nickel", "Penny", "Dime" };
    int coin_values[5] = { 100, 25, 5, 1, 10 };
    int cents = 0;
    for (size_t i = 0; i < 5; i++)
    {
        output << coins[i] << " " << coin_names[i] << endl;
        cout << coins[i] << " " << coin_names[i] << endl;
        cents += coin_values[i] * coins[i];
    }
    cout << "Total Amount: $" << cents / 100 << "." << (((cents % 100) >= 10) ? (to_string(cents % 100)) : (string("0").append(to_string(cents % 100)))) << endl;
    cout << cents << " Cents" << endl;
    output << "Total Amount: $" << cents / 100 << "." << (((cents % 100) >= 10) ? (to_string(cents % 100)) : (string("0").append(to_string(cents % 100)))) << endl;
}
int main(int argc, char** argv)//try doing A* pathing between maximums using distance + something related to the pixel intensity so that it tries to follow edges with missing pixels
{
    REVERSE_BINARIZATION = false;
    MAX_OFFSET = 20;
    MAX_CENTER_OFFSET = 8;
    THRESHOLD = 5000;
    THRESHOLD2 = 10000;
    
    //MAGNITUDE_THRESHOLD = 10;
    //N = 60;
    //part1();
    //part2();
    auto time = high_resolution_clock::now();
    part2(argc, argv);
    cout << std::chrono::duration_cast<seconds>(high_resolution_clock::now() - time).count();

    //line l(vec2d(0, 0), vec2d(1, 1));
    //cout << l.left_right(vec2d(0, 1)) << endl;
    //cout << l.left_right(vec2d(1, 0));
}

//TODO: Try using canny edge detection instead of image binarization, maybe remove the recursive part of it or otherwise tweak it to make it faster?