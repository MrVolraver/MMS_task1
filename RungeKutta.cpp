#include "RungeKutta.h"

#include <iostream>
#include <vector>

// Метод Рунге-Кутты 4-го порядка для системы ОДУ
void rungeKutta4(std::vector<double> (*f)(double, const std::vector<double>&), std::vector<double>& y, double t0, double t1, double h) {
    double t = t0;
    while (t < t1)
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

        t += h;

        for (int i=0; i<y.size(); i++)
            std::cout << y[i] << "  ";
        std::cout << std::endl;
    }
}
