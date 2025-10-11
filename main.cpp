#include <iostream>
#include <vector>
#include <windows.h>
#include <iostream>
#include <fstream>

#define EPS 1e-5
// коэффициенты
#define M 1
double Kabs[2] = {1.029, 2.481};
double Kel[2]  = {0.670, 0.649};
double K1[2]   = {0.118, 0.949};
double K2[2]   = {1.981, 2.757};
double K11[2]  = {0.201, 2.900};
double K12[2]  = {0.295, 0.466};
double K13[2]  = {8.998, 1.376};
double K21[2]  = {0.138, 0.531};
double K22[2]  = {0.347, 0.869};
double K23[2]  = {0.913, 2.814};
double K31[2]  = {0.687, 1.490};
double K32[2]  = {0.580, 1.666};
double K33[2]  = {2.343, 7.060};

// Система ОДУ
std::vector<double> System(double t, const std::vector<double>& y)
{
    std::vector<double> dy_dt;
    dy_dt.resize(y.size());

    dy_dt[0] = -Kabs[M] * y[0];
    dy_dt[1] = Kabs[M] * y[0] - (K1[M]+Kel[M]) * y[1] + K2[M] * y[2];
    dy_dt[2] = K1[M] * y[1] - K2[M] * y[2];
    dy_dt[3] = K11[M] * y[1] - K12[M] * y[3] + K13[M];
    dy_dt[4] = K21[M] * y[1] - K22[M] * y[4] + K23[M];
    dy_dt[5] = K31[M] * y[1] - K32[M] * y[5] + K33[M];

    return dy_dt;
}


// Метод Рунге-Кутты 4-го порядка для системы ОДУ
void rungeKutta4(std::vector<double> (*f)(double, const std::vector<double>&), std::vector<double>& y, double t0, double t1, double h) {
    double t = t0;
    int flag;
    while (t < t1)
    {
        flag = 1;
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
            if (delta > EPS)
                flag = 0;
        }

        if (flag) break;

        t += h;

        for (int i=0; i<y.size(); i++)
            std::cout << y[i] << "  ";
        std::cout << std::endl;
    }
}

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


int main()
{
    std::vector<double> y = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double t0 = 0.0, t1 = 10.0, h = 0.1;

    //Adams_Bashforth_Moulton(System, y, t0, t1, h);
    //rungeKutta4(System, y, t0, t1, h);
    Milne_Simpson(System, y, t0, t1, h);
    return 0;
}