// ConsoleApplication1.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <vector>
#include <fstream>

struct Point
{
    Point(int64_t X, int64_t Y) : x(X), y(Y) {}

    int64_t x = 0;
    int64_t y = 0;
};
 
int64_t defineDistance(Point & a, Point & b)
{
    int64_t dist = 0;

    return dist;
}

void getPointsFromFile(std::vector<Point>& points, size_t & amountPoints, std::ifstream & fin)
{
    for (size_t i = 0; i < amountPoints; ++i)
    {
        points.push_back(Point(0, 0));
        fin >> points[i].x >> points[i].y;
    }
}

int main()
{
    std::ifstream fin;
    std::ofstream fout;
    size_t amountPoints = 0;

    fin.open("input.txt");
    fin >> amountPoints;

    std::vector<std::vector<double>> graph(amountPoints, std::vector<double>(amountPoints));

    std::vector<Point> points;
    getPointsFromFile(points, amountPoints, fin);
    for (size_t i = 0; i < amountPoints; ++i)
    {
        std::cout << points[i].x << " " << points[i].y << "\n";
    }
    

    
    fout.open("output.txt");
    
    std::cout << "Hello World!\n";
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
