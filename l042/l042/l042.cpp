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
#include <deque>

using namespace std;
using namespace std::chrono;

int PPM_HEIGHT;
int PPM_WIDTH;
int N;
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
const pixel BLUE = { 0, 0, 255 };
inline void plot_pixel(int y, int x, pixel** pixels, pixel color = BLACK)
{
    //if (y >= 0 && x >= 0 && y < 400 && x < 400)
    if (((y | x) & 0x40000000) == 0 && y < 400 && x < 400)//bit trickery to reduce comparisons, I'm checking to see if either x or y have their negative bit set
        pixels[y][x] = color;
}
void draw_pixels(pixel** pixels)//writes pixels to a ppm file
{
    std::ofstream file("./grahamscan.ppm");
    if (!file.is_open())//in case something odd happens
    {
        cout << "FILE NOT OPENING";
        exit(1);
    }
    file.write("P3 400 400 255 \n", 16);

    for (size_t i = 0; i < 400; i++)
    {
        for (size_t j = 0; j < 400; j++)
        {
            file << std::to_string(pixels[i][j].r) << " " << std::to_string(pixels[i][j].g) << " " << std::to_string(pixels[i][j].b) << " ";
        }
        file << "\n";
    }
    file.flush();
    file.close();
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
};
void draw_horizontal_line(vec2i start, vec2i end, pixel** pixels)
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
    plot_pixel(end.get_y(), end.get_x(), pixels);
}
void draw_vertical_line(vec2i start, vec2i end, pixel** pixels)
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
    plot_pixel(end.get_y(), end.get_x(), pixels);
}
void draw_line(vec2i start, vec2i end, pixel** pixels)//fills in the end point specifically
{
    if (abs(end.get_y() - start.get_y()) <= abs(end.get_x() - start.get_x()))
        draw_horizontal_line(start, end, pixels);
    else
        draw_vertical_line(start, end, pixels);
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
void draw_incircle(vec2d p0, vec2d p1, vec2d p2, pixel** pixels)
{
    double sides[3];
    sides[0] = ((p1 - p2).absolute()).mag();//gets three side lengths, side[n] is opposite of points[n]
    sides[1] = ((p0 - p2).absolute()).mag();
    sides[2] = ((p0 - p1).absolute()).mag();
    double a = sides[0], b = sides[1], c = sides[2];//to reduce clutter
    double s = (a + b + c) / 2;
    vec2d center = { (a * p0.get_x() + b * p1.get_x() + c * p2.get_x()) / (a + b + c), (a * p0.get_y() + b * p1.get_y() + c * p2.get_y()) / (a + b + c) };//uses wikipedia algorithm based off of Heron's formula 
    double radius = round(sqrt((s - a) * (s - b) * (s - c) / s) * 400.);
    draw_circle_brehnsham(center * 400, radius, pixels);
}
vec2d draw_circumcircle(vec2d p0, vec2d p1, vec2d p2, pixel** pixels)
{
    vec2d intercept;//algorithm I got off wikipedia
    double D = 2 * (p0.get_x() * (p1.get_y() - p2.get_y()) + p1.get_x() * (p2.get_y() - p0.get_y()) + p2.get_x() * (p0.get_y() - p1.get_y()));
    intercept.set_x((1 / D) * (p0.mag2() * (p1.get_y() - p2.get_y()) + p1.mag2() * (p2.get_y() - p0.get_y()) + p2.mag2() * (p0.get_y() - p1.get_y())));
    intercept.set_y((1 / D) * (p0.mag2() * (p2.get_x() - p1.get_x()) + p1.mag2() * (p0.get_x() - p2.get_x()) + p2.mag2() * (p1.get_x() - p0.get_x())));
    double radius = round(intercept.dist(p0) * 400.);
    draw_circle_brehnsham(intercept * 400, radius, pixels);
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
    draw_line(vec2d(0, b), vec2d(400, slope1 * 400 + b), pixels);//could be more efficient but it doesn't really matter
}
pixel** pixels = new pixel * [400];

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
        draw_line(perp.p1 * 400, perp.p2 * 400, pixels);
        //draw_circle_brehnsham(p * 400, 7, pixels, BLUE);
        draw_line(p1 * 400, p2 * 400, pixels);
        return perp;
    }
    double distance_to_point(vec2d v)//returns squared value since I only use it for comparisons
    {
        return abs((p2.get_x() - p1.get_x()) * (p1.get_y() - v.get_y()) - (p1.get_x() - v.get_x()) * (p2.get_y() - p1.get_y())) / length;
    }
};

struct point_hasher
{
    size_t operator()(vec2d point) const
    {
        return (size_t)point.get_x() ^ (size_t)point.get_y();
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
//void part2()
//{
//    std::ofstream file("./points.txt");
//    file << setprecision(23);
//
//    ofstream results("./results.txt");
//    results << setprecision(23);
//    results << N << " points" << endl;
//    //results << "brute force" << endl;
//
//    vector<vec2d> points;
//    random_device dev;
//    std::mt19937_64 rng(dev());
//    std::uniform_real_distribution<> random_real(0, 1);
//    for (size_t i = 0; i < N; i++)
//    {
//        points.push_back({ random_real(rng), random_real(rng) });
//        file << points.back().get_x() << "  " << points.back().get_y() << endl;
//    }
//    for (size_t i = 0; i < 400; i++)
//    {
//        pixels[i] = new pixel[400];
//        memset(pixels[i], 255, sizeof(pixel) * 400);
//    }
//    for (vec2d v : points)
//        draw_circle_brehnsham(v * 400, 3, pixels);
//    //find min y coord
//    int min_y = 0;
//    for (size_t i = 0; i < N; i++)
//    {
//        if (points[i].get_y() < points[min_y].get_y())
//            min_y = i;
//    }
//    vec2d lowest_y = points[min_y];
//    points.erase(points.begin() + min_y);
//    sort(points.begin(), points.end(), [=](vec2d p, vec2d p2) {return line(lowest_y, p).get_angle() > line(lowest_y, p2).get_angle(); });
//
//    
//
//    deque<vec2d> convex = deque<vec2d>(points.begin(), points.end());
//    convex.push_front(lowest_y);
//    for (size_t i = 2; i < convex.size(); i++)
//    {
//        if (line(convex[i - 2], convex[i - 1]).left_right(convex[i]) > 0)//right turn so go back
//        {
//            i -= 1;
//            convex.erase(convex.begin() + i - 1);
//        }
//    }
//
//
//    for (vec2d v : convex)
//            draw_circle_brehnsham(v * 400, 5, pixels, RED);
//        //cout << v.get_x() << " : " << v.get_y() << endl;
//    //sort(convex.begin(), convex.end(), angle_pred);
//    for (size_t i = 0; i < convex.size() - 1; i++)
//    {
//        draw_line(convex[i] * 400, convex[i + 1] * 400, pixels);
//    }
//    draw_line(convex[0] * 400, convex[convex.size() - 1] * 400, pixels);
//    //draw_line(l.get_p1() * 400, l.get_p2() * 400, pixels);
//    draw_pixels(pixels);
//    cout << convex.size();
//}
//int main()
//{
//    N = 60;
//    part2();
//    //line l(vec2d(0, 0), vec2d(1, 1));
//    //cout << l.left_right(vec2d(0, 1)) << endl;
//    //cout << l.left_right(vec2d(1, 0));
//}