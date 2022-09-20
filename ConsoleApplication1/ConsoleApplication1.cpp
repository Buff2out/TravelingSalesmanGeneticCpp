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

void swapGens(size_t i, size_t j1, size_t j2, std::vector<std::vector<size_t>>& popul)
{
    size_t temp = popul[i][j1];
    popul[i][j1] = popul[i][j2];
    popul[i][j2] = temp;
}

void swapChromo(size_t i1, size_t i2, std::vector<std::vector<size_t>>& popul)
{
    std::vector<size_t> temp;
    for (size_t j = 0; j < popul[i1].size(); ++j)
    {
        temp.push_back(popul[i1][j]);
    }
    for (size_t j = 0; j < popul[i1].size(); ++j)
    {
        popul[i1][j] = popul[i2][j];
        popul[i2][j] = temp[j];
    }
}

void swapDists(size_t i1, size_t i2, std::vector<double>& dists)
{
    double temp = dists[i1];
    dists[i1] = dists[i2];
    dists[i2] = temp;
}

size_t partition(std::vector<double>& dists, std::vector<std::vector<size_t>>& popul, int64_t low, int64_t high, double pivot) {
    int64_t i = low;
    int64_t j = low;
    while (i <= high) 
    {
        if (dists[i] > pivot)
        {
            ++i;
        }
        else 
        {
            swapDists(i, j, dists);
            swapChromo(i, j, popul);
            ++i;
            ++j;
        }
    }
    return j - 1;
}

void quickSort(std::vector<double>& dists, std::vector<std::vector<size_t>>&popul, int64_t low, int64_t high) {
    if (low < high) 
    {
        double pivot = dists[high];
        int64_t pos = partition(dists, popul, low, high, pivot);
        quickSort(dists, popul, low, pos - 1);
        quickSort(dists, popul, pos + 1, high);
    }
}

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
        raNum = 1 + rand() % (amountPoints - 1); // -1 (длина хромосомы)
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

void toMutate(std::vector<std::vector<size_t>>& popul, size_t const& amountPoints, size_t const& k)
{
    // расслабон, по сравнению со скрещиванием просто инвертируем последовательность элементов в случайном сгенеринном диапазоне
    size_t a1 = 0, b1 = 0;
    for (size_t i = 0; i < 2 * (amountPoints - k); ++i)
    {
        a1 = rand() % (amountPoints - 2);
        b1 = a1 + rand() % (amountPoints - 1 - a1);
        for (size_t j = 0; j <= ((b1 - a1) / 2); ++j)
        {
            swapGens(i, a1 + j, b1 - j, popul);
        }
    }
}

double defineDist(std::vector<std::vector<double>>& graph, std::vector<size_t> chromoGen, size_t const& amountPoints, size_t const& k)
{
    double summ = 0;
    summ = summ + graph[0][chromoGen[0]];
    for (size_t j = 0; j < amountPoints - 2; ++j)
    {
        summ = summ + graph[chromoGen[j]][chromoGen[j + 1]];
    }
    summ = summ + graph[0][chromoGen[amountPoints - 2]];
    return summ;
}

void selectionAndSort(std::vector<std::vector<double>>& graph, std::vector<std::vector<size_t>>& popul, size_t const& amountPoints, size_t const& k)
{
    // а вот селекция это уже серьёзно
    std::vector<double> dists;
    for (size_t i = 0; i < 2 * (amountPoints - k); ++i)
    {
        dists.push_back(defineDist(graph, popul[i], amountPoints, k));
    }
    quickSort(dists, popul, 0, 2 * (amountPoints - k) - 1);
    // отбрасываем "невыживших" (а точнее этого даже делать не надо - всё потом само перепишется на моменте скрещивания
    // check selection
    for (size_t i = 0; i < 2 * (amountPoints - k); ++i)
    {
        std::cout << "Расстояние: " << dists[i] << ". Хромосома: ";
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

    getPointsFromFile(points, amountPoints, fin);
    
    fillEmptyMtrx(graph, amountPoints);
    fillDistsToMtrx(graph, points, amountPoints);

    // отсюда начинается зацикливание геналгоритма
    crossOver(popul, amountPoints, k);
    toMutate(popul, amountPoints, k);
    selectionAndSort(graph, popul, amountPoints, k);
    

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
