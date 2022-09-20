// ConsoleApplication1.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <iomanip>

struct Point
{
    Point(double X, double Y) : x(X), y(Y) {}

    double x = 0;
    double y = 0;
};

void swapGens(long int i, long int j1, long int j2, std::vector<std::vector<long int>>& popul)
{
    long int temp = popul[i][j1];
    popul[i][j1] = popul[i][j2];
    popul[i][j2] = temp;
}

void swapChromo(long int i1, long int i2, std::vector<std::vector<long int>>& popul)
{
    std::vector<long int> temp;
    for (long int j = 0; j < popul[i1].size(); ++j)
    {
        temp.push_back(popul[i1][j]);
    }
    for (long int j = 0; j < popul[i1].size(); ++j)
    {
        popul[i1][j] = popul[i2][j];
        popul[i2][j] = temp[j];
    }
}

void swapDists(long int i1, long int i2, std::vector<double>& dists)
{
    double temp = dists[i1];
    dists[i1] = dists[i2];
    dists[i2] = temp;
}

long int partition(std::vector<double>& dists, std::vector<std::vector<long int>>& popul, long int low, long int high, double pivot) {
    long int i = low;
    long int j = low;
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

void quickSort(std::vector<double>& dists, std::vector<std::vector<long int>>&popul, long int low, long int high) {
    if (low < high) 
    {
        double pivot = dists[high];
        long int pos = partition(dists, popul, low, high, pivot);
        quickSort(dists, popul, low, pos - 1);
        quickSort(dists, popul, pos + 1, high);
    }
}

void getPointsFromFile(std::vector<Point>& points, long int const & amountPoints, std::ifstream & fin)
{
    for (long int i = 0; i < amountPoints; ++i)
    {
        points.push_back(Point(0, 0));
        fin >> points[i].x >> points[i].y;
    }
}

void fillEmptyMtrx(std::vector<std::vector<double>>& graph, long int const & amountPoints)
{
    for (long int i = 0; i < amountPoints; ++i)
    {
        for (long int j = 0; j < amountPoints; ++j)
        {
            graph[i][j] = -1;
        }
    }
    for (long int j = 0; j < amountPoints; ++j)
    {
        graph[j][j] = 0;
    }
}

void fillDistsToMtrx(std::vector<std::vector<double>>& graph, std::vector<Point>& points, long int const& amountPoints)
{
    for (long int i = 0; i < amountPoints; ++i)
    {
        for (long int j = i + 1; j < amountPoints; ++j)
        {
            graph[j][i] = sqrt(pow((points[j].x - points[i].x), 2) + pow((points[j].y - points[i].y), 2));
            graph[i][j] = graph[j][i];
        }
    }
}

void createStartPop(std::vector<std::vector<long int>>& popul, long int const& amountPoints, long int const & k)
{
    for (long int i = 0; i < amountPoints - k; ++i)
    {
        for (long int j = 0; j < amountPoints - 1; ++j)
        {
            popul[i][j] = j + 1;
        }
    }
    for (long int i = 0; i < amountPoints; ++i)
    {
        std::random_shuffle(popul[i].begin(), popul[i].end());
    }
    //// // check popul
    //for (long int i = 0; i < amountPoints - k; ++i)
    //{
    //    for (long int j = 0; j < amountPoints - 1; ++j)
    //    {
    //        std::cout << popul[i][j] << " ";
    //    }
    //    std::cout << "\n";
    //}
}

void crossOver(std::vector<std::vector<long int>>& popul, long int const& amountPoints, long int const& k)
{
    std::vector<std::vector<long int>>::iterator popIt = popul.begin();
    //������� ������� ��������� ����� �������� ����������� �� �������� ������� popul
    for (long int i = 0; i < (amountPoints - k); ++i)
    {
        ++popIt;
    }
    // ������������ ������� ���������
    std::random_shuffle(popul.begin(), popIt);
    //���������� �������������� ��������� ����� ����� � ��������� ��������
    std::vector<bool> genBool1; // �� ����� ���� ��� ����������� ���? ��������� �� ���-�����
    std::vector<bool> genBool2; // �� ����� ���� ��� ����������� ���? ��������� �� ���-�����
    long int raNum = 0;
    // ��������� ������� ������ "��������"
    for (long int j = 0; j < amountPoints - 1; ++j) { genBool1.push_back(true); }
    for (long int j = 0; j < amountPoints - 1; ++j) { genBool2.push_back(true); }

    for (long int i = 0; i < amountPoints - k; i = i + 2)
    {
        raNum = 1 + rand() % (amountPoints - 1); // -1 (����� ���������)
        long int j1 = 0;
        // ������ raNum ����� ��������� � ��������
        while (j1 < raNum)
        {
            popul[amountPoints - k + i][j1] = popul[i][j1];
            genBool1[popul[i][j1] - 1] = false;
            popul[amountPoints - k + i + 1][j1] = popul[i + 1][j1];
            genBool2[popul[i + 1][j1] - 1] = false;
            ++j1;
        }
        // ��������� [raNum; amountPoints - 1) ����� ��������� � ���� �� ��������
        long int j2 = 0;
        long int j3 = j1;
        while (j1 < amountPoints - 1)
        {
            if (genBool1[popul[i + 1][j2] - 1])
            {
                popul[amountPoints - k + i][j1] = popul[i + 1][j2];
                ++j1;
            }
            ++j2;
        }
        // ���������� ����� ����, �� �� ������ ��������
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
        for (long int j = 0; j < amountPoints - 1; ++j) { genBool1[j] = true; }
        for (long int j = 0; j < amountPoints - 1; ++j) { genBool2[j] = true; }
    }


    //// // check popul
    //for (long int i = 0; i < 2 * (amountPoints - k); ++i)
    //{
    //    for (long int j = 0; j < amountPoints - 1; ++j)
    //    {
    //        std::cout << popul[i][j] << " ";
    //    }
    //    std::cout << "\n";
    //}
}

void toMutate(std::vector<std::vector<long int>>& popul, long int const& amountPoints, long int const& k)
{
    // ���������, �� ��������� �� ������������ ������ ����������� ������������������ ��������� � ��������� ����������� ���������
    long int a1 = 0, b1 = 0;
    for (long int i = 0; i < 2 * (amountPoints - k); ++i)
    {
        if (amountPoints - 2 == 0)
        {
            a1 = 0;
        }
        else
        {
            a1 = rand() % (amountPoints - 2);
        }
        b1 = a1 + rand() % (amountPoints - 1 - a1);
        for (long int j = 0; j <= ((b1 - a1) / 2); ++j)
        {
            swapGens(i, a1 + j, b1 - j, popul);
        }
    }
}

double defineDist(std::vector<std::vector<double>>& graph, std::vector<long int> chromoGen, long int const& amountPoints, long int const& k)
{
    double summ = 0;
    summ = summ + graph[0][chromoGen[0]];
    for (long int j = 0; j < amountPoints - 2; ++j)
    {
        summ = summ + graph[chromoGen[j]][chromoGen[j + 1]];
    }
    summ = summ + graph[0][chromoGen[amountPoints - 2]];
    return summ;
}

void selectionAndSort(std::vector<double>& dists, std::vector<std::vector<double>>& graph, std::vector<std::vector<long int>>& popul, long int const& amountPoints, long int const& k)
{
    // � ��� �������� ��� ��� ��������
    for (long int i = 0; i < 2 * (amountPoints - k); ++i)
    {
        dists[i] = defineDist(graph, popul[i], amountPoints, k);
    }
    quickSort(dists, popul, 0, 2 * (amountPoints - k) - 1);
    // ����������� "����������" (� ������ ����� ���� ������ �� ���� - �� ����� ���� ����������� �� ������� �����������
    //// check selection
    //for (long int i = 0; i < 2 * (amountPoints - k); ++i)
    //{
    //    std::cout << "����������: " << dists[i] << ". ���������: ";
    //    for (long int j = 0; j < amountPoints - 1; ++j)
    //    {
    //        std::cout << popul[i][j] << " ";
    //    }
    //    std::cout << "\n";
    //}
}

int main()
{
    std::ifstream fin;
    std::ofstream fout;

    long int amount = 0;
    fin.open("input.txt");
    fin >> amount;
    long int const amountPoints = amount;

    std::vector<std::vector<double>> graph(amountPoints, std::vector<double>(amountPoints));

    long int k = 0;
    if (amountPoints % 2 == 1)
    {
        k = 1;
    }
    std::vector<std::vector<long int>> popul(2 * (amountPoints - k), std::vector<long int>(amountPoints - 1));
    std::vector<Point> points;
    std::vector<double> dists;
    for (long int i = 0; i < 2 * (amountPoints - k); ++i)
    {
        dists.push_back(0);
    }

    createStartPop(popul, amountPoints, k);

    getPointsFromFile(points, amountPoints, fin);

    fillEmptyMtrx(graph, amountPoints);
    fillDistsToMtrx(graph, points, amountPoints);

    long int summ = 0;
    double theBest = -1;
    while (summ != 3)
    {
        crossOver(popul, amountPoints, k);
        toMutate(popul, amountPoints, k);
        selectionAndSort(dists, graph, popul, amountPoints, k);
        if (theBest == dists[0])
        {
            ++summ;
        }
        else
        {
            theBest = dists[0];
            summ = 0;
        }
    }
    fout.open("output.txt");
    fout << std::setprecision(6) << std::fixed << dists[0] << "\n";
    fout << 1 << " ";
    for (long int j = 0; j < amountPoints - 1; ++j)
    {
        fout << popul[0][j] + 1 << " ";
    }
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
