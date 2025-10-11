#include "Milne_Simpson.h"

#include <iostream>
#include <vector>

// Метод Милна-Симсона
void Milne_Simpson(std::vector<double> (*f)(double, const std::vector<double>&), std::vector<double>& y, double t0, double t1, double h)
{
    double t = t0;
    std::vector<double> yim3;
    yim3.resize(y.size());
    std::vector<double> yim2;
    yim2.resize(y.size());
    std::vector<double> yim1;
    yim1.resize(y.size());
    std::vector<double> yi;
    yi.resize(y.size());

    yim3=y;

    // Рунге-Кутта для первых 3х шагов
    for (int i=0; i<3; i++)
    {
        std::vector<double> k1 = f(t, y);
        std::vector<double> yk1;
        for (int i=0; i<y.size(); i++)
            yk1.push_back(y[i]+h*k1[i]/2);
        std::vector<double> k2 = f(t + h / 2, yk1);
        std::vector<double> yk2;
        for (int i=0; i<y.size(); i++)
            yk2.push_back(y[i]+h*k2[i]/2);
        std::vector<double> k3 = f(t + h / 2, yk2);
        std::vector<double> yk3;
        for (int i=0; i<y.size(); i++)
            yk3.push_back(y[i]+h*k2[i]/2);
        std::vector<double> k4 = f(t + h, yk3);

        for (int i=0; i<y.size(); i++)
        {
            double delta =  h / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
            y[i] += delta;
        }

        switch (i)
        {
        case 0:
            {
                yim2=y;
                break;
            }
        case 1:
            {
                yim1=y;
                break;
            }
        }
        t += h;
    }
    yi = y;

    // Милн-Симпсон
    while (t < t1)
    {
        //предиктор - Милн
        std::vector<double> ki = f(t, yi);
        std::vector<double> kiM1 = f(t-h, yim1);
        std::vector<double> kiM2 = f(t-2*h, yim2);
        std::vector<double> y_;
        y_.resize(y.size());

        for (int i=0; i<y.size(); i++)
        {
            y_[i] = yim3[i] + 4*h/3 * (2*ki[i] - kiM1[i] + 2*kiM2[i]);
        }

        std::vector<double> kiP1 = f(t+h, y_);

        for (int i=0; i<y.size(); i++)
        {
            y[i] = yim1[i] + h/3 * (kiP1[i] + 4*ki[i] + kiM1[i]);
        }
        t+=h;

        yim3 = yim2;
        yim2 = yim1;
        yim1 = yi;
        yi = y;

        for (int i=0; i<y.size(); i++)
            std::cout << y[i] << "  ";
        std::cout << std::endl;

    }
}