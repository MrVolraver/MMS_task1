#include "Adams_Bashforth_Moulton.h"

#include <iostream>
#include <vector>

// Метод Адамса-Башфорта-Моултона порядка для системы ОДУ
void Adams_Bashforth_Moulton(std::vector<double> (*f)(double, const std::vector<double>&), std::vector<double>& y, double t0, double t1, double h) {
    double t = t0;
    std::vector<double> y_1;
    y_1.resize(y.size());
    std::vector<double> y_2;
    y_2.resize(y.size());
    std::vector<double> y_3;
    y_3.resize(y.size());

    y_1=y;
 // Первые 4 шага вычисляются при помощи Рунге-Кутты
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
                y_2=y;
                break;
            }
        case 1:
            {
                y_3=y;
                break;
            }
        }
        t += h;
    }

    while (t < t1)
    {
        //предиктор - Башфорт
        std::vector<double> k1 = f(t, y);
        std::vector<double> k2 = f(t-h, y_3);
        std::vector<double> k3 = f(t-2*h, y_2);
        std::vector<double> k4 = f(t - 3*h, y_1);
        std::vector<double> y_;
        y_.resize(y.size());
        for (int i=0; i<y.size(); i++)
        {
            y_[i] = y[i] + h/24 * (55*k1[i] - 59*k2[i] + 37*k3[i] - 9*k4[i]);
        }

        k1 = f(t+h, y_);
        k2 = f(t, y);
        k3 = f(t - h, y_3);
        k4 = f(t - 2*h, y_2);

        y_1 = y_2;
        y_2 = y_3;
        y_3 = y;
        // Корректор - Молутон
        for (int i=0; i<y.size(); i++)
        {
            y[i] += h/24 * (9*k1[i] + 19*k2[i] - 5*k3[i] + k4[i]);
        }
        t+=h;

        for (int i=0; i<y.size(); i++)
            std::cout << y[i] << "  ";
        std::cout << std::endl;

    }
}