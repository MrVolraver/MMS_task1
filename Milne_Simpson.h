#pragma once
#include <vector>

extern "C" __declspec(dllexport) void Milne_Simpson(std::vector<double> (*f)(double, const std::vector<double>&), std::vector<double>&, double, double, double);