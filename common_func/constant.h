#pragma once
#include <complex>
#include <numbers>

using namespace std;

// int N = 100;
// int steps = 4000 * 1.05;
// int M=2*N+1;
// double dx = 2 * numbers::pi / M;
complex<double> i_unit(0.0, 1.0);

double L = 2 * numbers::pi;
int N = 256;
int M = N;
int steps = 100;
double dx = L / N;
double alpha = 0.01;
double dt = 1.00 / steps;