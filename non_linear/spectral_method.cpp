#include <iostream>
#include <numbers>
#include <Eigen/Dense>
#include <fstream>
#include <algorithm>
#include <complex>
#include "fourier.h"
#include "constant.h"
#include "runge_kutta.h"

using namespace std;

Eigen::VectorXcd non_linear_term(Eigen::VectorXcd &u_hat)
{
  Eigen::VectorXcd sum = Eigen::VectorXcd::Zero(M);
  for (int m = 0; m < M; m++)
  {
    int k_wave = m - N;
    for (int l = 0; l < M; l++)
    {
      int p = l - N;
      int q = k_wave - p;
      int m_p = q + N;
      if (m_p >= 0 && m_p < M)
        sum(m) -= i_unit * static_cast<double>(q) * u_hat(l) * u_hat(m_p);
    }
  }

  return sum;
}

int main()
{
  ofstream file("non_linear/spectral_non_linear.csv");
  Eigen::VectorXcd u(M), du1(M), du2(M), du3(M), du4(M), k1(M), k2(M), k3(M), k4(M), tmp(M);
  u.setZero();
  for (int i = 0; i < M; i++)
    u[i] = sin(i * dx);

  Eigen::VectorXcd u_hat = forward_fourier(u);

  for (int i = 0; i < steps; i++)
  {
    u_hat = runge_kutta(u_hat, non_linear_term);
  }
  u = inverse_fourier(u_hat);
  for (int i = 0; i < M; i++)
  {
    file << i * dx << "," << u(i).real() << "\n";
  }

  cout << "計算が完了しました！" << endl;
  return 0;
}