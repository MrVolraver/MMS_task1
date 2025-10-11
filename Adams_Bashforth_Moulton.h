#pragma once
#include <vector>

extern "C" __declspec(dllexport) void Adams_Bashforth_Moulton(std::vector<double> (*f)(double, const std::vector<double>&), std::vector<double>&, double, double, double);