#pragma once
#include <vector>
#include <iostream>

extern "C" __declspec(dllexport) void RungeKutta(std::vector<double>(double, std::vector<double>&), std::vector<double>&, double, double, double, std::ofstream*);