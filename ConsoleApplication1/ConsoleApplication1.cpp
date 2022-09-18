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

void createStartPop(std::vector<std::vector<size_t>>& popul, size_t const& amountPoints, size_t const & k)
{
    for (size_t i = 0; i < amountPoints - k; ++i)
    {
        for (size_t j = 0; j < amountPoints - 1; ++j)
        {
            popul[i][j] = j + 1;
        }
    }
    for (size_t i = 0; i < amountPoints; ++i)
    {
        std::random_shuffle(popul[i].begin(), popul[i].end());
    }
    //// // check popul
    //for (size_t i = 0; i < amountPoints - k; ++i)
    //{
    //    for (size_t j = 0; j < amountPoints - 1; ++j)
    //    {
    //        std::cout << popul[i][j] << " ";
    //    }
    //    std::cout << "\n";
    //}
}

void crossOver(std::vector<std::vector<size_t>>& popul, size_t const& amountPoints, size_t const& k)
{
    std::vector<std::vector<size_t>>::iterator popIt = popul.begin();
    //костыль который позволяет найти итератор указывающий на середину вектора popul
    for (size_t i = 0; i < (amountPoints - k); ++i)
    {
        ++popIt;
    }
    // перемешиваем будущих родителей
    std::random_shuffle(popul.begin(), popIt);
    //скрещиваем перемешавшихся родителей между собой и порождаем потомков
    std::vector<bool> genBool1; // на самом деле это хромосомный бул? состоящий из бул-генов
    std::vector<bool> genBool2; // на самом деле это хромосомный бул? состоящий из бул-генов
    size_t raNum = 0;
    // заполняем булевый список "трушками"
    for (size_t j = 0; j < amountPoints - 1; ++j) { genBool1.push_back(true); }
    for (size_t j = 0; j < amountPoints - 1; ++j) { genBool2.push_back(true); }

    for (size_t i = 0; i < amountPoints - k; i = i + 2)
    {
        raNum = 1 + rand() % (amountPoints - 1); // -1 (длина хромосомы) -1 (последний элемент)
        size_t j1 = 0;
        // первые raNum генов добавляем в потомков
        while (j1 < raNum)
        {
            popul[amountPoints - k + i][j1] = popul[i][j1];
            genBool1[popul[i][j1] - 1] = false;
            popul[amountPoints - k + i + 1][j1] = popul[i + 1][j1];
            genBool2[popul[i + 1][j1] - 1] = false;
            ++j1;
        }
        // следующие [raNum; amountPoints - 1) генов добавляем в этих же потомков
        size_t j2 = 0;
        size_t j3 = j1;
        while (j1 < amountPoints - 1)
        {
            if (genBool1[popul[i + 1][j2] - 1])
            {
                popul[amountPoints - k + i][j1] = popul[i + 1][j2];
                ++j1;
            }
            ++j2;
        }
        // повторение всего того, но со вторым потомком
        j2 = 0;
        j1 = j3;
        while (j1 < amountPoints - 1)
        {
            if (genBool2[popul[i][j2] - 1])
            {
                popul[amountPoints - k + i + 1][j1] = popul[i][j2];
                ++j1;
            }
            ++j2;
        }
        for (size_t j = 0; j < amountPoints - 1; ++j) { genBool1[j] = true; }
        for (size_t j = 0; j < amountPoints - 1; ++j) { genBool2[j] = true; }
    }


    // // check popul
    for (size_t i = 0; i < 2 * (amountPoints - k); ++i)
    {
        for (size_t j = 0; j < amountPoints - 1; ++j)
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

    size_t k = 0;
    if (amountPoints % 2 == 1)
    {
        k = 1;
    }
    std::vector<std::vector<size_t>> popul(2 * (amountPoints - k), std::vector<size_t>(amountPoints - 1));
    std::vector<Point> points;

    createStartPop(popul, amountPoints, k);
    crossOver(popul, amountPoints, k);
    getPointsFromFile(points, amountPoints, fin);
    
    fillEmptyMtrx(graph, amountPoints);
    fillDistsToMtrx(graph, points, amountPoints);
    //std::vector<int64_t> answerSeq(amountPoints);


    //// // check graph
    //for (size_t i = 0; i < amountPoints; ++i)
    //{
    //    for (size_t j = 0; j < amountPoints; ++j)
    //    {
    //        std::cout << graph[i][j] << " ";
    //    }
    //    std::cout << "\n";
    //}
    
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
