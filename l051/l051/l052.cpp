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

using namespace std;
using namespace std::chrono;

int IMG_HEIGHT;
int IMG_WIDTH;
int N;
int THRESHOLD2;
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
    if (((y | x) & 0x40000000) == 0 && y < IMG_HEIGHT && x < IMG_WIDTH)//bit trickery to reduce comparisons, I'm checking to see if either x or y have their negative bit set
        pixels[y][x] = color;
}
void draw_pixels(pixel** pixels, string name)//writes pixels to a ppm file
{
    std::ofstream file("./" + name);
    if (!file.is_open())//in case something odd happens
    {
        cout << "FILE NOT OPENING";
        exit(1);
    }
    file.write(("P3 " + to_string(IMG_WIDTH) + " " + to_string(IMG_HEIGHT) + " 255 \n").c_str(), 10 + to_string(IMG_HEIGHT).size() + to_string(IMG_WIDTH).size());

    for (size_t i = 0; i < IMG_HEIGHT; i++)
    {
        for (size_t j = 0; j < IMG_WIDTH; j++)
        {
            file << std::to_string(pixels[i][j].r) << " " << std::to_string(pixels[i][j].g) << " " << std::to_string(pixels[i][j].b) << " ";
        }
        file << "\n";
    }
}
void draw_pixels(int* pixels, string name)//writes pixels to a ppm file
{
    std::ofstream file("./" + name);
    if (!file.is_open())//in case something odd happens
    {
        cout << "FILE NOT OPENING";
        exit(1);
    }
    file.write(("P3 " + to_string(IMG_WIDTH) + " " + to_string(IMG_HEIGHT) + " 255 \n").c_str(), 10 + to_string(IMG_HEIGHT).size() + to_string(IMG_WIDTH).size());

    for (size_t i = 0; i < IMG_HEIGHT * IMG_WIDTH; i++)
    {
        file << std::to_string(pixels[i]) << " " << std::to_string(pixels[i]) << " " << std::to_string(pixels[i]) << " ";
    }
    file << "\n";
    
}
void draw_pixels(int** pixels, string name)//writes pixels to a ppm file
{
    std::ofstream file("./" + name);
    if (!file.is_open())//in case something odd happens
    {
        cout << "FILE NOT OPENING";
        exit(1);
    }
    file.write(("P3 " + to_string(IMG_WIDTH) + " " + to_string(IMG_HEIGHT) + " 255 \n").c_str(), 10 + to_string(IMG_HEIGHT).size() + to_string(IMG_WIDTH).size());

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
        return (22801763489 + (point.get_x())) * 22801763489 + (point.get_y());//magic number is bilionth prime number
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
        int mat [3][3];
        int row1 [3];
        int row2 [3];
        int row3 [3];
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
    void kernel(int16_t* pixels, int16_t** kernel)// evaluate a kernel on an image
    {
        for (int y = 2; y < IMG_HEIGHT - 1; y++)//ignore edges and sobel will be initialized black
        {
            for (int x = 2; x < IMG_WIDTH - 1; x++)
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
                kernel[y][x] = val;
            }
        }

    }

};

void part1(bool gaussian = false)//maybe try remaking with a color sobel filter, one per channel
{
    ifstream ppm("image.ppm");
    string line;
    getline(ppm, line);
    getline(ppm, line);
    IMG_WIDTH = stoi(line.substr(0, line.find(' ')));
    IMG_HEIGHT = stoi(line.substr(line.find(' ') + 1));
    getline(ppm, line);

    int* pixels = new int [IMG_HEIGHT * IMG_WIDTH];//black and white
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
        while(pos < line.size())
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
    draw_pixels(pixels, "imageg.ppm");
    //const 
    //int* pixels_gaus = new int[IMG_HEIGHT * IMG_WIDTH];//black and white gausianed if requested
    //for (int y = 1; y < IMG_HEIGHT - 1; y++)
    //{
    //    for (int x = 1; x < IMG_WIDTH - 1; x++)
    //    {
    //        pixels_gaus[y * IMG_WIDTH + IMG_HEIGHT] = ()
    //    }
    //}

    int** sobelh = new int * [IMG_HEIGHT];
    for (size_t i = 0; i < IMG_HEIGHT; i++)
    {
        sobelh[i] = new int[IMG_WIDTH];
        memset(sobelh[i], 0, IMG_WIDTH * sizeof(int));
    }
    int** sobelv = new int * [IMG_HEIGHT];
    for (size_t i = 0; i < IMG_HEIGHT; i++)
    {
        sobelv[i] = new int[IMG_WIDTH];
        memset(sobelv[i], 0, IMG_WIDTH * sizeof(int));
    }
    int** sobel = new int * [IMG_HEIGHT];
    for (size_t i = 0; i < IMG_HEIGHT; i++)
    {
        sobel[i] = new int[IMG_WIDTH];
        memset(sobel[i], 0, IMG_WIDTH * sizeof(int));
    }

    int SOBELH[3][3] = { {1, 0, -1}, {2, 0, -2}, {1, 0, -1} };
    int SOBELV[3][3] = { {1, 2, 1}, {0, 0, 0}, {-1, -2, -1} };
    mat3x3 mat_sobel_horiz = mat3x3(SOBELH);
    mat3x3 mat_sobel_vert = mat3x3(SOBELV);

    mat_sobel_horiz.kernel(pixels, sobelh);
    mat_sobel_vert.kernel(pixels, sobelv);

    for (int y = 0; y < IMG_HEIGHT; y++)
    {
        for (int x = 0; x < IMG_WIDTH; x++)
        {
            //int val = ;
            sobel[y][x] = round(sqrt(sobelh[y][x] * sobelh[y][x] + sobelv[y][x] * sobelv[y][x]));//min(val > THRESHOLD ? val : 0, 255);
        }
    }
    draw_pixels(sobel, "imagem.ppm");
}
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
            if((((pos.get_y() + y > 0) && (pos.get_x() + x > 0))) && (y < IMG_HEIGHT) && (x < IMG_WIDTH) && (thresh[pos.get_y()+y][pos.get_x()+x] > 0))
                flood_fill(vec2i(pos.get_x() + x, pos.get_y() + y), thresh, depth + 1);
        }
    }
}
void part2()
{
    ifstream ppm("image.ppm");
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
    draw_pixels(pixels, "imageg.ppm");
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

    int SOBELH[3][3] = { {1, 0, -1}, {2, 0, -2}, {1, 0, -1} };
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

    double scale = 255. / (double)max;
    for (int y = 0; y < IMG_HEIGHT; y++)
    {
        for (int x = 0; x < IMG_WIDTH; x++)
        {
            sobel[y][x] = min(255., sobel[y][x] * scale);//min(val > THRESHOLD ? val : 0, 255);
        }
    }
    draw_pixels(sobel, "imagem.ppm");

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
        }
    }
    draw_pixels(sobel, "image1.ppm");
}
int main()//try doing A* pathing between maximums using distance + something related to the pixel intensity so that it tries to follow edges with missing pixels
{
    THRESHOLD = 20;
    THRESHOLD2 = 100;//lower than might be expected because I dynamically scale the combined Sobel values to make the maximum as close to 255 as possible resulting in significantly less values attempting to go above 255
    //N = 60;
    //part1();
    part2();
    //line l(vec2d(0, 0), vec2d(1, 1));
    //cout << l.left_right(vec2d(0, 1)) << endl;
    //cout << l.left_right(vec2d(1, 0));
}