// ConsoleApplication1.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <algorithm>

struct Point
{
    Point(double X, double Y) : x(X), y(Y) {}

    double x = 0;
    double y = 0;
};

void getPointsFromFile(std::vector<Point>& points, size_t const & amountPoints, std::ifstream & fin)
{
    for (size_t i = 0; i < amountPoints; ++i)
    {
        points.push_back(Point(0, 0));
        fin >> points[i].x >> points[i].y;
    }
}

void fillEmptyMtrx(std::vector<std::vector<double>>& graph, size_t const & amountPoints)
{
    for (size_t i = 0; i < amountPoints; ++i)
    {
        for (size_t j = 0; j < amountPoints; ++j)
        {
            graph[i][j] = -1;
        }
    }
    for (size_t j = 0; j < amountPoints; ++j)
    {
        graph[j][j] = 0;
    }
}

void fillDistsToMtrx(std::vector<std::vector<double>>& graph, std::vector<Point>& points, size_t const& amountPoints)
{
    for (size_t i = 0; i < amountPoints; ++i)
    {
        for (size_t j = i + 1; j < amountPoints; ++j)
        {
            graph[j][i] = sqrt(pow((points[j].x - points[i].x), 2) + pow((points[j].y - points[i].y), 2));
            graph[i][j] = graph[j][i];
        }
    }
}

void createStartPop(std::vector<std::vector<size_t>>& popul, size_t const& amountPoints)
{
    for (size_t i = 0; i < amountPoints; ++i)
    {
        for (size_t j = 0; j < amountPoints; ++j)
        {
            popul[i][j] = j;
        }
    }
    for (size_t i = 0; i < amountPoints; ++i)
    {
        std::random_shuffle(popul[i].begin(), popul[i].end());
    }
    // // check
    for (size_t i = 0; i < amountPoints; ++i)
    {
        for (size_t j = 0; j < amountPoints; ++j)
        {
            std::cout << popul[i][j] << " ";
        }
        std::cout << "\n";
    }
}

int main()
{
    std::ifstream fin;
    std::ofstream fout;

    size_t amount = 0;
    fin.open("input.txt");
    fin >> amount;
    size_t const amountPoints = amount;

    std::vector<std::vector<double>> graph(amountPoints, std::vector<double>(amountPoints));
    std::vector<std::vector<size_t>> popul(amountPoints, std::vector<size_t>(amountPoints));

    createStartPop(popul, amountPoints);
    //std::vector<int64_t> answerSeq(amountPoints);

    std::vector<Point> points;
    getPointsFromFile(points, amountPoints, fin);
    
    fillEmptyMtrx(graph, amountPoints);
    fillDistsToMtrx(graph, points, amountPoints);



    for (size_t i = 0; i < amountPoints; ++i)
    {
        for (size_t j = 0; j < amountPoints; ++j)
        {
            std::cout << graph[i][j] << " ";
        }
        std::cout << "\n";
    }
    
    //fout.open("output.txt");
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
