//#include "stdafx.h"
#include <fstream>
#include <cstdio>
#include <cstddef>
#include <ctime>
#include <iostream>
#include <array>
#include <vector>
#include <cmath>
#include <string>
#include <list>
#include <float.h>
#include <iomanip>
#include <ctime>
#include <algorithm>
#include <chrono>
#include <unordered_map>
#include <stack>
#include <unordered_set>

using namespace std;
int numPoints = 60;
int dim = 400;
int width, height;
double* point;

class Pixel {
private:
    int r = 0;
    int g = 0;
    int b = 0;
public:
    void setRGB(int r, int g, int b) {
        this->r = r;
        this->g = g;
        this->b = b;
    }

    void illuminate(string color) {
        if (color == "Black") {
            this->r = 0;
            this->g = 0;
            this->b = 0;
        }
        else if (color == "Red") {
            this->r = 255;
            this->g = 0;
            this->b = 0;
        }
        else if (color == "White") {
            this->r = 255;
            this->g = 255;
            this->b = 255;
        }
    }

    int getR() { return r; }
    int getG() { return g; }
    int getB() { return b; }
    int getRGBSum() { return r + g + b; }
    void print() {
        cout << r << " " << g << " " << b << " ";
    }
};

class PPMImage {
public:
    string PPMType; //P3 or P6
    int width;
    int height;
    int maxColorVal;
    Pixel** board;

    PPMImage() {
        ;
    }

    Pixel** getBoard() {
        return board;
    }
};

class Point {
    int x = 0;
    int y = 0;
public:
    Point() {}

    Point(int i, int j) {
        x = i; y = j;
    }

    void print() {
        cout << x << ' ' << y << endl;
    }

    int getX() { return x; }
    int getY() { return y; }

    bool operator == (const Point& P) const {
        return (x == P.x) && (y == P.y);
    }
};

void print(vector<Point> v) {
    for (size_t j = 0; j < v.size(); j++) {
        v[j].print();
        cout << " ";
    }
    cout << endl;
}

PPMImage readPPM(std::ifstream& in, PPMImage& img) {
    in >> img.PPMType;

    string w;
    in >> w;
    width = atoi(w.data());
    img.width = width;

    string h;
    in >> h;
    height = atoi(h.data());
    img.height = height;

    string maxVal;
    in >> maxVal;
    img.maxColorVal = atoi(maxVal.data());

    img.board = new Pixel * [img.width];
    for (int i = 0; i < img.width; i++) {
        img.board[i] = new Pixel[img.height];
    }

    for (int j = 0; j < img.height; j++) {
        for (int i = 0; i < img.width; i++) {
            if (i == 0 || j == 0 || i == img.width - 1 || j == img.height - 1) {
                int r; int g; int b;
                in >> r; in >> g; in >> b;
                img.board[i][j].setRGB(0, 0, 0);
            }
            else {
                int r; int g; int b;
                in >> r; in >> g; in >> b;
                img.board[i][j].setRGB(r, g, b);
            }
        }
    }
    return img;
}

PPMImage greyscale() {
    ifstream in;
    in.open("grumpycat.ppm", std::ifstream::binary);
    PPMImage p;
    readPPM(in, p);
    in.close();

    //Greying
    for (int j = 0; j < p.height; j++) {
        for (int i = 0; i < p.width; i++) {
            int greyVal = round(p.board[i][j].getRGBSum() / 3);
            p.board[i][j].setRGB(greyVal, greyVal, greyVal);
        }
    }

    Pixel** board = p.getBoard();

    ofstream ofs("imageg.ppm", ios_base::out | ios_base::binary);
    ofs << "P3" << endl << p.width << ' ' << p.height << endl << "255" << endl;
    for (int j = 0; j < p.height; j++)
    {
        for (int i = 0; i < p.width; i++)
        {
            Pixel pix = board[i][j];
            ofs << to_string(pix.getR()) << " " << to_string(pix.getG()) << " " << to_string(pix.getB()) << " ";
        }
    }
    return p;
}

vector<Point> getTwos(Pixel** board) {
    vector<Point> v;
    for (int i = 1; i < width - 1; i++) {
        for (int j = 1; j < height - 1; j++) {
            //board[i][j].print();
            if (board[i][j].getR() == 2) {
                Point p(i, j);
                v.push_back(p);
            }
        }
    }
    return v;
}

bool inRange(Point P) {
    if (P.getX() == 0 || P.getX() == width - 1 || P.getY() == 0 || P.getY() == height - 1)
        return false;
    return true;
}

int int_hash(int key) {
    key = key ^ (key << 2);
    return key;
}

struct PointHash {
    int operator()(Point pnt) const {
        return (53 + pnt.getX()) * 53 + pnt.getY();
    }
};

unordered_set<Point, PointHash> seen;

void floodFill(Pixel** board, Point P, int depth) {
    if (depth > 1000) return;
    if (inRange(P) && seen.find(P) != seen.end()) {
        int x = P.getX();
        int y = P.getY();
        if (board[x][y].getR() == 1 || board[x][y].getR() == 2) {
            seen.insert(P);
            Point south(x, y - 1);
            floodFill(board, south, depth + 1);
            Point north(x, y + 1);
            floodFill(board, north, depth + 1);
            Point west(x - 1, y);
            floodFill(board, west, depth + 1);
            Point east(x + 1, y);
            floodFill(board, east, depth + 1);
            Point northeast(x + 1, y + 1);
            floodFill(board, northeast, depth + 1);
            Point northwest(x - 1, y + 1);
            floodFill(board, northwest, depth + 1);
            Point southeast(x + 1, y - 1);
            floodFill(board, southeast, depth + 1);
            Point southwest(x - 1, y - 1);
            floodFill(board, southwest, depth + 1);
        }
    }
    return;
}

void sobel(PPMImage img) {
    Pixel** board = img.board;

    Pixel** newBoard = new Pixel * [img.width];
    for (int i = 0; i < img.width; i++) {
        newBoard[i] = new Pixel[img.height];
    }

    for (int i = 1; i < img.width - 2; i++) {
        for (int j = 1; j < img.height - 2; j++) {
            int horiVal = (board[i][j].getR() + (2 * board[i][j + 1].getR()) + board[i][j + 2].getR()) -
                (board[i + 2][j].getR() + (2 * board[i + 2][j + 1].getR()) + board[i + 2][j + 2].getR());
            int vertVal = (board[i][j + 2].getR() + (2 * board[i + 1][j + 2].getR()) + board[i + 2][j + 2].getR()) -
                (board[i][j].getR() + (2 * board[i + 1][j].getR()) + board[i + 2][j].getR());
            int mag = abs(horiVal) + abs(vertVal);
            if (mag > 150)
                newBoard[i][j].setRGB(2, 2, 2);
            else if (mag > 70)
                newBoard[i][j].setRGB(1, 1, 1);
            else
                newBoard[i][j].setRGB(0, 0, 0);
        }
    }

    ofstream ofs1("imagesplit.ppm", ios_base::out | ios_base::binary);
    ofs1 << "P3" << endl << img.width << ' ' << img.height << endl << "2" << endl;
    for (int j = 0; j < img.height; j++)
    {
        for (int i = 0; i < img.width; i++)
        {
            Pixel pix = newBoard[i][j];
            ofs1 << to_string(pix.getR()) << " " << to_string(pix.getG()) << " " << to_string(pix.getB()) << " ";
        }
    }

    vector<Point> v = getTwos(newBoard); //Gets all twos

    for (Point p : v) { //For each two
        floodFill(newBoard, p, 0); //Fill the twos and connected 1's
    }

    for (int i = 0; i < width; i++) { //Nuke the board to black
        for (int j = 0; j < height; j++) {
            newBoard[i][j].setRGB(0, 0, 0);
        }
    }

    for (Point p : seen) { //Mark the seen 1's and 2's as 1's
        newBoard[p.getX()][p.getY()].setRGB(1, 1, 1);
    }

    ofstream ofs("imagem.ppm", ios_base::out | ios_base::binary);
    ofs << "P3" << endl << img.width << ' ' << img.height << endl << "1" << endl;
    for (int j = 0; j < img.height; j++)
    {
        for (int i = 0; i < img.width; i++)
        {
            Pixel pix = newBoard[i][j];
            ofs << to_string(pix.getR()) << " " << to_string(pix.getG()) << " " << to_string(pix.getB()) << " ";
        }
    }
}

void part1() {
    PPMImage img = greyscale();
    sobel(img);
}

void part2() {
    PPMImage img = greyscale();
    sobel(img);

}

int main() {
    //part1();
    part2();
}