#pragma once
#include <vector>
#include <iostream>

extern "C" __declspec(dllexport) void Adams_Bashforth_Moulton(std::vector<double>(double, std::vector<double>&), std::vector<double>&, double, double, double, std::ofstream*);