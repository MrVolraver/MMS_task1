#include <iostream>
#include <vector>
#include <windows.h>
#include <iostream>
#include <fstream>
#include <functional>

#define EPS 1e-6
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
std::vector<double> System(double t, std::vector<double>& y)
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

// метод Гаусса
std::vector<double> gauss(std::vector<std::vector<double>>& A, std::vector<double>& b) {
    size_t row_size = A.size();
    size_t col_size = A.back().size();
    // Прямой ход Гаусса
    double pivot = 0.0;
    for (size_t i = 0; i < row_size; i++) {
        for (size_t j = i + 1; j < col_size; j++) {
            if (std::abs(A.at(j).at(i)) < EPS) {
                continue;
            }
            pivot = A.at(j).at(i) / A.at(i).at(i);
            b.at(j) -= pivot * b.at(i);
            for (size_t k = 0; k < row_size; k++) {
                A.at(j).at(k) -= pivot * A.at(i).at(k);
            }
        }
    }

    // Обратный ход Гаусса
    std::vector<double> x(row_size);
    for (int i = row_size - 1; i >= 0; i--) {
        x.at(i) = b.at(i);
        for (size_t j = i + 1; j < row_size; j++) {
            x.at(i) -= x.at(j) * A.at(i).at(j);
        }
        if (A.at(i).at(i) == 0)
            return {-1,-1,-1,-1};
        x.at(i) /= A.at(i).at(i);
    }

    return x;
}

// аналитическое решение системы
std::vector<double> analytical_solution(double t, std::vector<std::vector<double>> initial)
{
    // Initial - начальные условия;   задаются как initial = {t, c}, где yi(t) = c
    // Решение системы
    std::vector<double> y;
    y.resize(6);

    // Решение первой подсистемы из y1 y2 y3 методом Эйлера
    // Собственные числа лямбда |A-labmdaE| = 0
    double lambda_1 = Kabs[M];
    double lambda_2 = ((-1)*(K1[M]+Kel[M]+K2[M]) + (sqrt((K1[M]+Kel[M]+K2[M])*(K1[M]+Kel[M]+K2[M])-4*K2[M]*Kel[M])))/2;
    double lambda_3 = ((-1)*(K1[M]+Kel[M]+K2[M]) - (sqrt((K1[M]+Kel[M]+K2[M])*(K1[M]+Kel[M]+K2[M])-4*K2[M]*Kel[M])))/2;

    //Собственные векторы для лямбд |A-labmdaE| = 0
    std::vector<double> e1;
    std::vector<double> e2;
    std::vector<double> e3;

    std::vector<std::vector<double>> A1 = {
        {-Kabs[M]-lambda_1,         0,                           0},
        {Kabs[M],                   -(K1[M]+Kel[M])-lambda_1,    K2[M]},
        {0,                         K1[M],                      -K2[M]-lambda_1}
    };
    std::vector<std::vector<double>> A2 = {
        {-Kabs[M]-lambda_2,  0,                  0},
        {Kabs[M],   -(K1[M]+Kel[M])-lambda_2,    K2[M]},
        {0,         K1[M],              -K2[M]-lambda_2}
    };
    std::vector<std::vector<double>> A3 = {
        {-Kabs[M]-lambda_3,  0,                  0},
        {Kabs[M],   -(K1[M]+Kel[M])-lambda_3,    K2[M]},
        {0,         K1[M],              -K2[M]-lambda_3}
    };

    std::vector<double> b = {0,0,0};

    e1 = gauss(A1,b);
    e2 = gauss(A3,b);
    e3 = gauss(A3,b);

    // вычисление С из начальных условий RC = I для у
    std::vector<double> C;

    std::vector<std::vector<double>> R = {
        {exp(lambda_1*initial[0][0])*e1[0],  exp(lambda_2*initial[0][0])*e2[0], exp(lambda_3*initial[0][0])*e3[0]},
        {exp(lambda_1*initial[1][0])*e1[1],exp(lambda_2*initial[1][0])*e2[1],exp(lambda_3*initial[1][0])*e3[1]},
        {exp(lambda_1*initial[2][0])*e1[2],exp(lambda_2*initial[2][0])*e2[2],exp(lambda_3*initial[2][0])*e3[2]}
    };

    std::vector<double> I = {initial[0][1], initial[1][1],initial[2][1]};

    C = gauss(R,I);

    for (int i=0; i<3;i++)
        y[i] = C[0]*exp(lambda_1*t)*e1[i] + C[1]*exp(lambda_2*t)*e2[i] + C[2]*exp(lambda_3*t)*e3[i];

    // вычисление С из начальных условий RC = I для x
    double C4 = (initial[3][1] - (K13[M]-K11[M]*y[1])/(K12[M]) * exp(K12[M]*initial[3][0]));
    double C5 = (initial[4][1] - (K23[M]-K21[M]*y[1])/(K22[M]) * exp(K22[M]*initial[4][0]));
    double C6 = (initial[5][1] - (K33[M]-K31[M]*y[1])/(K32[M]) * exp(K32[M]*initial[5][0]));

    y[3] = C4/exp(K12[M]*t) + (K13[M]-K11[M]*y[1])/(K12[M]);
    y[4] = C5/exp(K22[M]*t) + (K23[M]-K21[M]*y[1])/(K22[M]);
    y[5] = C6/exp(K32[M]*t) + (K33[M]-K31[M]*y[1])/(K32[M]);

    return y;

}

int main()
{
    std::vector<double> y = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double t0 = 0.0, t1 = 10.0, h = 0.1;

    HMODULE dll1 = LoadLibrary("RungeKutta.dll");
    typedef int FT1(std::vector<double>(double, std::vector<double>&), std::vector<double>&, double, double, double, std::ofstream*);
    FT1* RungeKutta = (FT1*)GetProcAddress(dll1, "RungeKutta");

    HMODULE dll2 = LoadLibrary("Adams_Bashforth_Moulton.dll");
    typedef int FT2(std::vector<double>(double, std::vector<double>&), std::vector<double>&, double, double, double, std::ofstream*);
    FT2* Adams_Bashforth_Moulton = (FT2*)GetProcAddress(dll2, "Adams_Bashforth_Moulton");

    HMODULE dll3 = LoadLibrary("Milne_Simpson.dll");
    typedef int FT3(std::vector<double>(double, std::vector<double>&), std::vector<double>&, double, double, double, std::ofstream*);
    FT3* Milne_Simpson = (FT3*)GetProcAddress(dll3, "Milne_Simpson");

    char* name = (char*)"RungeKutta.csv";
    std::ofstream csv1;
    csv1.open(name);

    char* name2 = (char*)"Milne_Simpson.csv";
    std::ofstream csv2;
    csv2.open(name2);

    char* name3 = (char*)"Adams_Bashforth_Moulton.csv";
    std::ofstream csv3;
    csv3.open(name3);

    RungeKutta(System, y, t0, t1, h, &csv1);
    y = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    Adams_Bashforth_Moulton(System, y, t0, t1, h, &csv3);
    y = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    Milne_Simpson(System, y, t0, t1, h, &csv2);

    // начальные условия
    std::vector<std::vector<double>> initial = {
            {0,0},
            {0,0},
            {0,0},
            {0,0},
            {0,0},
            {0,0},
    };

    char* name4 = (char*)"Analytical_solution.csv";
    std::ofstream csv4;
    csv4.open(name4);

    for (double t=0.0; t<t1+h; t+=h)
    {
        std::vector<double> ya = analytical_solution(t, initial);
        for (int i=0; i<5; i++)
            csv4 << ya[i] << ", " ;
        csv4 << ya[5] << std::endl;
    }
    return 0;
}