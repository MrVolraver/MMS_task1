#include <iostream>
#include <vector>
#include <windows.h>
#include <iostream>
#include <fstream>
#include <functional>

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

void func(double a)
{
    std::cout << K31[M] << std::endl;
}

int main()
{
    std::vector<double> y = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double t0 = 0.0, t1 = 10.0, h = 0.1;

    HMODULE dll1 = LoadLibrary("RungeKutta.dll");
    typedef int FT1(std::vector<double>(double, std::vector<double>&), std::vector<double>&, double, double, double);
    FT1* RungeKutta = (FT1*)GetProcAddress(dll1, "RungeKutta");

    HMODULE dll2 = LoadLibrary("Adams_Bashforth_Moulton.dll");
    typedef int FT2(std::vector<double>(double, std::vector<double>&), std::vector<double>&, double, double, double);
    FT2* Adams_Bashforth_Moulton = (FT2*)GetProcAddress(dll2, "Adams_Bashforth_Moulton");

    HMODULE dll3 = LoadLibrary("Milne_Simpson.dll");
    typedef int FT3(std::vector<double>(double, std::vector<double>&), std::vector<double>&, double, double, double);
    FT3* Milne_Simpson = (FT3*)GetProcAddress(dll3, "Milne_Simpson");

    RungeKutta(System, y, t0, t1, h);
    //Adams_Bashforth_Moulton(System, y, t0, t1, h);
    //Milne_Simpson(System, y, t0, t1, h);

    return 0;
}