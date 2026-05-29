#pragma once
#include <complex>

using namespace std;

int N = 100;
int M = 2 * N + 1;
int steps = 4000;
double dt = 1.05 / steps;
double dx = 2 * M_PI / M;
complex<double> i_unit(0.0, 1.0);