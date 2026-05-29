#pragma once
#include <Eigen/Dense>
#include <complex>
#include "constant.h"

using namespace std;

Eigen::VectorXcd forward_fourier(Eigen::VectorXcd &u)
{
    Eigen::VectorXcd u_hat = Eigen::VectorXcd::Zero(M);
    for (int i = 0; i < M; i++)
    {
        int k_wave = i - N;
        for (int j = 0; j < M; j++)
        {
            double x = j * dx;
            u_hat(i) += u(j) * exp(-i_unit * static_cast<double>(k_wave) * x);
        }
        u_hat(i) /= M;
    }
    return u_hat;
}

Eigen::VectorXcd inverse_fourier(Eigen::VectorXcd &u_hat)
{
    Eigen::VectorXcd u = Eigen::VectorXcd::Zero(M);
    for (int i = 0; i < M; i++)
    {
        double x = i * dx;
        for (int j = 0; j < M; j++)
        {
            int k_wave = j - N;
            u(i) += u_hat(j) * exp(i_unit * static_cast<double>(k_wave) * x);
        }
    }
    return u;
}
