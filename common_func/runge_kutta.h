#pragma once
#include <Eigen/Dense>
#include "constant.h"

template <typename VecType, typename F>
VecType runge_kutta(VecType &u, F f)
{
    VecType du1(M), du2(M), du3(M), du4(M), k1(M), k2(M), k3(M), k4(M), tmp(M);

    du1 = f(u);
    k1 = dt * du1;

    tmp = u + k1 / 2;

    du2 = f(tmp);
    k2 = dt * du2;

    tmp = u + k2 / 2;

    du3 = f(tmp);
    k3 = dt * du3;

    tmp = u + k3;

    du4 = f(tmp);
    k4 = dt * du4;

    u += (k1 + 2 * k2 + 2 * k3 + k4) / 6;
    return u;
}