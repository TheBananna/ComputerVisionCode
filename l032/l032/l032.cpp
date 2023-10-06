//Nicholas Bonanno PD: 3
#define _USE_MATH_DEFINES
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
#include <vector>
#include <list>
#include <algorithm>
#include <tuple>

using namespace std;
using namespace std::chrono;


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
    bool operator <(const vec2d& v2)
    {
        return x < v2.x;
    }
    void print()
    {
        cout << x << " : " << y << endl;
    }
};
vec2d* get_points()
{
    vec2d* points = new vec2d[4];
    for (size_t i = 0; i < 4; i++)
    {
        points[i] = { (double)rand() / RAND_MAX, (double)rand() / RAND_MAX };
    }
    return points;
}
struct pixel//pixel, I'm using byte to reduce memory usage. This brings a measureable speedup, ~10% on my laptop
{
public:
    uint8_t r, g, b;
};
const pixel BLACK = { 0, 0, 0 };
const pixel WHITE = { 255, 255, 255 };
const pixel RED = { 255, 0, 0 };
inline void plot_pixel(int y, int x, pixel** pixels, pixel color = BLACK)
{
    //if (y >= 0 && x >= 0 && y < 800 && x < 800)
    if (((y | x) & 0x80000000) == 0 && y < 800 && x < 800)//bit trickery to reduce comparisons, I'm checking to see if either x or y have their negative bit set
        pixels[y][x] = color;
}
void draw_pixels(pixel** pixels)//writes pixels to a ppm file
{
    std::ofstream file("./points.ppm");
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
    double radius = round(sqrt((s - a) * (s - b) * (s - c) / s) * 800.);
    draw_circle_brehnsham(center * 800, radius, pixels);
}
vec2d draw_circumcircle(vec2d p0, vec2d p1, vec2d p2, pixel** pixels)
{
    vec2d intercept;//algorithm I got off wikipedia
    double D = 2 * (p0.get_x() * (p1.get_y() - p2.get_y()) + p1.get_x() * (p2.get_y() - p0.get_y()) + p2.get_x() * (p0.get_y() - p1.get_y()));
    intercept.set_x((1 / D) * (p0.mag2() * (p1.get_y() - p2.get_y()) + p1.mag2() * (p2.get_y() - p0.get_y()) + p2.mag2() * (p0.get_y() - p1.get_y())));
    intercept.set_y((1 / D) * (p0.mag2() * (p2.get_x() - p1.get_x()) + p1.mag2() * (p0.get_x() - p2.get_x()) + p2.mag2() * (p1.get_x() - p0.get_x())));
    double radius = round(intercept.dist(p0) * 800.);
    draw_circle_brehnsham(intercept * 800, radius, pixels);
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
    double slope1 = p1.slope(p2);//uses y = mx + b to allow easy calculation of start and end points, any pixels drawn outside of the 800x800 square will be discarded
    double b = p1.get_y() - slope1 * p1.get_x();
    draw_line(vec2d(0, b), vec2d(800, slope1 * 800 + b), pixels);//could be more efficient but it doesn't really matter
}
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
        draw_line_to_edge(p1 * 800, p2 * 800, pixels);
        //draw_circle_brehnsham(p1 * 800, 5, pixels);
        //draw_circle_brehnsham(p2 * 800, 5, pixels);

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
    line get_perpendicular(vec2d p)
    {
        line perp;
        perp.angle = angle + M_PI_2 * (p.get_y() < y_at(p.get_x()) ? 1 : -1);
        perp.length = length;
        perp.p1 = p;
        perp.p2 = p + (vec2d(cos(perp.angle), sin(perp.angle)) * length);

        vec2d intersect = intersection(perp);
        perp.p1 = intersect;
        perp.p2 = intersect + (vec2d(cos(perp.angle), sin(perp.angle)) * intersect.dist(p1));
        perp.length = intersect.dist(perp.p2);
        return perp;
    }
};
vec2d get_nth(list<vec2d> l, int n)
{
    int i = 0;
    for (vec2d& v : l)
    {
        if (i == n)
            return v;
        i++;
    }
    return vec2d();
}
vec2d* brute_force(vector<vec2d>& ps, int N = 60)
{
    list<vec2d> points;
    for (size_t i = 0; i < N; i++)
    {
        points.push_back(ps[i]);
    }
    vec2d* closest = new vec2d[2];
    closest[0] = get_nth(points, 0);
    closest[1] = get_nth(points, 1);
    double dist = closest[0].dist2(closest[1]);

    int adv = 0;
    for (vec2d& p : points)
    {
        adv++;
        auto p2 = points.begin();
        advance(p2, adv);
        for (; p2 != points.end(); ++p2)
        {
            double dist_new = p.dist2(*p2);
            if (dist > dist_new)
            {
                closest[0] = p;
                closest[1] = *p2;
                dist = dist_new;
            }
        }
    }
    return closest;
}
void part1()
{
    //cout << "start" << endl;
    //high_resolution_clock::time_point t1 = high_resolution_clock::now();
    const int N = 60;

    std::ofstream file("./points.txt");
    file << setprecision(23);


    list<vec2d> points;
    for (size_t i = 0; i < N; i++)
    {
        points.push_back(vec2d::rand_vec());
        file << points.back().get_x() << "  " << points.back().get_y() << endl;
    }
    vec2d* closest = new vec2d[2];
    closest[0] = get_nth(points, 0);
    closest[1] = get_nth(points, 1);
    double dist = closest[0].dist2(closest[1]);

    int adv = 0;
    for (vec2d& p : points)
    {
        adv++;
        auto p2 = points.begin();
        advance(p2, adv);
        for (; p2 != points.end(); ++p2)
        {
            double dist_new = p.dist2(*p2);
            if (dist > dist_new)
            {
                closest[0] = p;
                closest[1] = *p2;
                dist = dist_new;
            }
        }
    }

    //high_resolution_clock::time_point t2 = high_resolution_clock::now();
    //duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
    //std::cout << time_span.count() << " seconds" << endl;

    cout << dist << endl;
    cout << closest[0].get_x() << "  " << closest[0].get_y() << endl;
    cout << closest[1].get_x() << "  " << closest[1].get_y() << endl;
    pixel** pixels = new pixel * [800];//inits pixel grid to white
    for (size_t i = 0; i < 800; i++)
    {
        pixels[i] = new pixel[800];
        memset(pixels[i], 255, 800 * sizeof(pixel));
    }
    for (vec2d& point : points)
    {
        draw_circle_brehnsham(point * 800, 3, pixels);
    }

    draw_circle_brehnsham(closest[0] * 800, 3, pixels, RED);
    draw_circle_brehnsham(closest[1] * 800, 3, pixels, RED);


    draw_pixels(pixels);
}
char* scratch = NULL;
char* current_scratch;
char* stop_scratch;
constexpr uint32_t SCRATCH_SIZE = 1024 * 1024 * 128;//64 MiB
void init_scratch()
{
    if (scratch)
        delete[] scratch;
    scratch = new char[SCRATCH_SIZE];
    current_scratch = scratch;
    stop_scratch = scratch + SCRATCH_SIZE;
}

void* get_scratch(char size)
{
    if (current_scratch + size > stop_scratch)
        current_scratch = stop_scratch - SCRATCH_SIZE;
    char* to_ret = current_scratch;
    current_scratch += size;
    return to_ret;
}
long points_in_range(vector<vec2d>* points, int mid, int min, int max, double d)
{
    vec2d p_middle = points->at(mid);
    int new_min = min;
    int new_max = max;
    while (p_middle.get_x() - d > points->at(new_min).get_x())
        new_min++;
    while (p_middle.get_x() + d < points->at(new_max).get_x())
        new_max--;
    return new_min | (new_max >> 32);
}
int* strip_left(vector<vec2d>* points, int mid, int min, double d)
{
    int pos = mid - 1;
    vec2d midval = points->at(mid);
    while (pos >= min && midval.dist2(points->at(pos)) < d)//binary search possibility?
        pos--;
    pos++;
    int* ret_points = (int*)get_scratch(sizeof(int) * (mid - pos + 1));
    ret_points[0] = mid - pos;
    for (int i = pos; i < mid; i++)
    {
        ret_points[i + 1] = i;
    }
    return ret_points;
}
int* strip_right(vector<vec2d>* points, int mid, int max, double d)
{
    int pos = mid + 1;
    vec2d midval = points->at(mid);
    while (pos <= max && midval.dist2(points->at(pos)) < d)//binary search possibility?
        pos++;
    pos--;
    int* ret_points = (int*)get_scratch(sizeof(int) * (pos - mid + 1));
    ret_points[0] = pos - mid;
    for (int i = mid + 1; i < pos + 1; i++)
    {
        ret_points[i + 1] = i;
    }
    return ret_points;
}
std::tuple<double, int, int> zipper_merge(vector<vec2d>* points, int min, int max)
{
    if (max - min == 2)
    {
        return std::tuple<double, int, int> (points->at(max).dist2(points->at(min)), min, max);
    }
    if (max - min == 3)
    {
        double dist1 = points->at(min).dist2(points->at(min + 1));
        double dist2 = points->at(min).dist2(points->at(max));
        double dist3 = points->at(min + 1).dist2(points->at(max));
        if (dist1 < dist2)
        {
            if (dist1 < dist3)
            {
                return std::make_tuple(dist1, min, min+1);
            }
            return std::make_tuple(dist3, min+1, max);
        }
        if (dist2 < dist3)
            return std::make_tuple(dist2, min, max);
        return std::make_tuple(dist3, min + 1, max);
    }

    auto left = zipper_merge(points, min, (max + min) / 2);
    auto right = zipper_merge(points, (max + min) / 2 + 1, max);
    auto smaller = get<0>(left) < get<0>(right) ? left : right;
    double d = get<0>(smaller);

    int middle = (min + max) / 2;
    int* left_strip = strip_left(points, middle, min, d);
    int* right_strip = strip_right(points, middle, max, d);

    for (size_t i = 0; i < left_strip[0]; i++)
    {
        vec2d left_p = points->at(i);
        for (size_t j = 0; j < right_strip[0]; j++)
        {
            double dist = left_p.dist2(points->at(j));//maybe make temp of points->at(j) so that two memory lookups aren't needed
            if (dist < d)
            {
                d = dist;
                smaller = std::make_tuple(dist, i, j);
            }
        }
    }

    return smaller;
}

void part2()
{
    //cout << "start" << endl;
    
    const int N = 60;
    init_scratch();
    std::ifstream file("./points.txt");
    string s;

    vector<vec2d> points;
    for (size_t i = 0; i < N; i++)
    {
        std::getline(file, s);
        int pos = s.find(' ');
        points.push_back(vec2d(stod(s.substr(0, pos)), stod(s.substr(pos + 2, s.length() - pos - 2))));
        //file << points.back().get_x() << "  " << points.back().get_y() << endl;
        //points.push_back(vec2d::rand_vec());
    }
    //cout << "sort start" << endl;
    //high_resolution_clock::time_point t1 = high_resolution_clock::now();
    std::sort(points.begin(), points.end());

    auto close = zipper_merge(&points, 0, points.size() - 1);
    cout << (sqrt(get<0>(close))) << endl;
    points[get<1>(close)].print();
    points[get<2>(close)].print();


    //cout << "Smallest point ( " << points[0].get_x() << " , " << points[0].get_y() << " )" << endl;

    //high_resolution_clock::time_point t2 = high_resolution_clock::now();
    //duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
    //std::cout << time_span.count() << " seconds" << endl;

    
    //pixel** pixels = new pixel * [800];//inits pixel grid to white
    //for (size_t i = 0; i < 800; i++)
    //{
    //    pixels[i] = new pixel[800];
    //    memset(pixels[i], 255, 800 * sizeof(pixel));
    //}
    //for (vec2d& point : points)
    //{
    //    draw_circle_brehnsham(point * 800, 3, pixels);
    //}

    //draw_circle_brehnsham(closest[0] * 800, 3, pixels, RED);
    //draw_circle_brehnsham(closest[1] * 800, 3, pixels, RED);


    //draw_pixels(pixels);
}
int main()
{
    srand(time(0));
    srand(0);
    part1();
    part2();
}