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
inline void plot_pixel(int y, int x, pixel** pixels, pixel color = BLACK)
{
    //if (y >= 0 && x >= 0 && y < 800 && x < 800)
    if (((y | x) & 0x80000000) == 0 && y < 800 && x < 800)//bit trickery to reduce comparisons, I'm checking to see if either x or y have their negative bit set
        pixels[y][x] = color;
}
void draw_pixels(pixel** pixels)//writes pixels to a ppm file
{
    std::ofstream file("./output.ppm");
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
void draw_circle_brehnsham(vec2i center, int radius, pixel** pixels)
{
    plot_pixel(center.get_y(), center.get_x(), pixels);
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
        plot_pixel(center.get_y() + y, center.get_x() + x, pixels);
        plot_pixel(center.get_y() - y, center.get_x() + x, pixels);
        plot_pixel(center.get_y() + y, center.get_x() - x, pixels);
        plot_pixel(center.get_y() - y, center.get_x() - x, pixels);

        plot_pixel(center.get_y() + x, center.get_x() + y, pixels);
        plot_pixel(center.get_y() - x, center.get_x() + y, pixels);
        plot_pixel(center.get_y() + x, center.get_x() - y, pixels);
        plot_pixel(center.get_y() - x, center.get_x() - y, pixels);

    }
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
vec2d part1()
{
    vec2d* points = get_points();
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
        file << "(" << points[i].get_x() << "," << points[i].get_y() << ") , ";
    }
    file << "(" << point4.get_x() << "," << point4.get_y() << ")";
    file.flush();
    file.close();
    return point4;
}
class square
{
public:
    vec2d* points;
    double area;
    square min(square s)
    {
        if (area < s.area)
        {
            return *this;
        }
        else
        {
            return s;
        }
    }
};
square get_square_area(void** inputs)
{
    vec2d A = *(vec2d*)inputs[0], B = *(vec2d*)inputs[1], C = *(vec2d*)inputs[2], D = *(vec2d*)inputs[3];
    line AB = line(A, B);
    line perp_AB = AB.get_perpendicular(C);

    line perp_perp_AB = perp_AB.get_perpendicular(D);
    //if (reverse)
    //    perp_perp_AB = line(perp_perp_AB.get_p2(), perp_perp_AB.get_p1());
    line l3 = AB.get_perpendicular(perp_perp_AB.get_p2());

    /*AB.draw(pixels);
    perp_AB.draw(pixels);
    perp_perp_AB.draw(pixels);
    l3.draw(pixels);*/

    double side = AB.intersection(perp_AB).dist(l3.intersection(AB));
    square s;
    s.points = new vec2d[4];
    s.points[0] = A;
    s.points[1] = B;
    s.points[2] = C;
    s.points[3] = D;
    s.area = side * side;
    return s;
}
void print_square(vec2d* inputs, pixel** pixels)
{
    vec2d A = inputs[0], B = inputs[1], C = inputs[2], D = inputs[3];
    line AB = line(A, B);
    line perp_AB = AB.get_perpendicular(C);

    line perp_perp_AB = perp_AB.get_perpendicular(D);
    //if (reverse)
    //    perp_perp_AB = line(perp_perp_AB.get_p2(), perp_perp_AB.get_p1());
    line l3 = AB.get_perpendicular(perp_perp_AB.get_p2());

    AB.draw(pixels);
    perp_AB.draw(pixels);
    perp_perp_AB.draw(pixels);
    l3.draw(pixels);
}
void*** prep_inputs(vec2d* A, vec2d* B, vec2d* C, vec2d* D)
{
    void*** inputs = new void** [6];
    for (size_t i = 0; i < 6; i++)
    {
        for (size_t j = 0; j < 5; j++)
        {
            inputs[i] = new void*[5];
        }
    }
    inputs[0][0] = A;
    inputs[0][1] = B;
    inputs[0][2] = C;
    inputs[0][3] = D;
    inputs[0][4] = 0;

    inputs[1][0] = A;
    inputs[1][1] = C;
    inputs[1][2] = B;
    inputs[1][3] = D;
    inputs[1][4] = 0;

    inputs[2][0] = A;
    inputs[2][1] = D;
    inputs[2][2] = B;
    inputs[2][3] = C;
    inputs[2][4] = 0;

    inputs[3][0] = A;
    inputs[3][1] = B;
    inputs[3][2] = D;
    inputs[3][3] = C;
    inputs[3][4] = (void*)1;

    inputs[4][0] = A;
    inputs[4][1] = C;
    inputs[4][2] = D;
    inputs[4][3] = B;
    inputs[4][4] = (void*)1;

    inputs[5][0] = A;
    inputs[5][1] = D;
    inputs[5][2] = C;
    inputs[5][3] = B;
    inputs[5][4] = (void*)1;
    return inputs;
}
void part2()
{
    pixel** pixels = new pixel * [800];//inits pixel grid to white
    for (size_t i = 0; i < 800; i++)
    {
        pixels[i] = new pixel[800];
        memset(pixels[i], 255, 800 * sizeof(pixel));
    }

    //vec2d* points = get_points();
    //vec2d p0, p1, p2, p3;
    //p0 = points[0];
    //p1 = points[1];
    //p2 = points[2];
    //p3 = points[3];


    fstream f("points.txt");
    string s;
    string newstring;
    std::getline(f, s);
    cout << s << endl;
    for (size_t i = 0; i < s.length(); i++)
    {
        if ((s[i] != '(') && (s[i] != ')') && (s[i] != ' ') && (s[i] != ','))
            newstring += s[i];
    }
    s = newstring;
    string point_strings[8];
    double el_doubles[8];

    int start = 2;
    for (size_t i = 0; i < 7; i++)
    {
        int temp = start;
        while (s[temp] != '.')
            temp++;
        el_doubles[i] = stod(s.substr(start - 1, temp - (start - 1)));
        start = temp + 1;
    }
    el_doubles[7] = stod(s.substr(start - 1, s.length() - (start - 1)));
    int x = 0;
    vec2d* points = new vec2d[4];
    for (size_t i = 0; i < 4; i++)
    {
        points[i] = { el_doubles[2 * i], el_doubles[2 * i + 1] };
        //cout << points[i].get_x() << " : " << points[i].get_y() << endl;
    }
    vec2d A = points[0], B = points[1], C = points[2], D = points[3];

    void*** inputs = prep_inputs(&A, &B, &C, &D);



    //print_square(inputs[3], pixels);
    //draw_circle_brehnsham(perp_perp_AB.get_p1() * 800, 5, pixels);
    //draw_circle_brehnsham(perp_perp_AB.get_p2() * 800, 5, pixels);
    //draw_circle_brehnsham(intersect2 * 800, 5, pixels);
    for (size_t i = 0; i < 4; i++)
        draw_circle_brehnsham(points[i] * 800, 2, pixels);

    std::ofstream file("./output.txt");
    file << setprecision(18);
    for (size_t i = 0; i < 3; i++)
    {
        points[i] = { el_doubles[2 * i], el_doubles[2 * i + 1] };
        file << '(' << points[i].get_x() << "," << points[i].get_y() << ") , ";
    }
    file << '(' << points[3].get_x() << "," << points[3].get_y() << ")" << endl;

    square smallest = { NULL, 1000 };
    for (size_t i = 0; i < 6; i++)
    {
        square sq = get_square_area(inputs[i]);
        smallest = sq.min(smallest);
        cout << "(" << sq.points[0].get_x() << "," << sq.points[0].get_y() << ") , (" << sq.points[1].get_x() << "," << sq.points[1].get_y() << ") , (" << sq.points[2].get_x() << "," << sq.points[2].get_y() << ") , (" << sq.points[3].get_x() << "," << sq.points[3].get_y() << ") Area=" << sq.area << endl;
        file << "(" << sq.points[0].get_x() << "," << sq.points[0].get_y() << ") , (" << sq.points[1].get_x() << "," << sq.points[1].get_y() << ") , (" << sq.points[2].get_x() << "," << sq.points[2].get_y() << ") , (" << sq.points[3].get_x() << "," << sq.points[3].get_y() << ") Area=" << sq.area << endl;
    }

    print_square(smallest.points, pixels);
    if (!file.is_open())//in case something odd happens
    {
        cout << "FILE NOT OPENING";
        exit(1);
    }
    file.flush();
    file.close();
    draw_pixels(pixels);

}
int main()
{
    //high_resolution_clock::time_point t1 = high_resolution_clock::now();
    srand(time(0));
    //part1();
    part2();
    //print_square(inputs[sq.count], pixels);
    //for (size_t i = 0; i < 3; i++)
    //{
    //    draw_circle_brehnsham({ points[i].x * 800, points[i].y * 800 }, 5, pixels);
    //}
    //draw_circle_brehnsham({ point4.x * 800, point4.y * 800 }, 10, pixels);

    //draw_line(points[1] * 800, points[2] * 800, pixels);
    //draw_line(points[2] * 800, points[0] * 800, pixels);
    //draw_line(points[0] * 800, points[1] * 800, pixels);

    //high_resolution_clock::time_point t2 = high_resolution_clock::now();
    //duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
    //std::cout << time_span.count();


}

