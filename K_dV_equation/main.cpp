#include <numbers>
#include <Eigen/Dense>
#include <fftw3.h>
#include <complex>
#include <conversion_fftw_eigen.h>
#include <runge_kutta.h>
#include <fstream>
#include <constant.h>
#include <iostream>

using namespace std;

bool is_linear_deleted = true;

Eigen::VectorXcd conversion_method(Eigen::VectorXcd &u_hat)
{
  Eigen::VectorXcd G_hat(N), u(N), sum(N);
  sum.setZero();
  eigen_to_fftw(u_hat, in);
  fftw_execute(p_backward);

  for (int i = 0; i < N; i++)
  {
    u(i) = complex<double>(out[i][0], out[i][1]) / static_cast<double>(N);
  }
  for (int i = 0; i < N; i++)
  {
    in[i][0] = u(i).real() * u(i).real() / 2.0;
    in[i][1] = 0.0;
  }
  fftw_execute(p_forward);

  G_hat = fftw_to_eigen(out);

  Eigen::VectorXcd F_hat(N);
  for (int i = 0; i < N; i++)
  {
    double k = i - N / 2;
    F_hat(i) = i_unit * k * G_hat(i);
    sum(i) = -F_hat(i) + alpha * i_unit * k * k * k * u_hat(i);
  }
  for (int i = 0; i < N; i++)
  {
    int k_wave = i - N / 2;
    if (abs(k_wave) > N / 4)
      sum(i) = 0.0;
  }
  return sum;
}

Eigen::VectorXcd conversion_method_linear_delete(double t, Eigen::VectorXcd &u_hat)
{
  Eigen::VectorXcd G_hat(N), u(N), sum(N);
  sum.setZero();
  eigen_to_fftw(u_hat, in);
  fftw_execute(p_backward);

  for (int i = 0; i < N; i++)
  {
    u(i) = complex<double>(out[i][0], out[i][1]) / static_cast<double>(N);
  }
  for (int i = 0; i < N; i++)
  {
    in[i][0] = u(i).real() * u(i).real() / 2.0;
    in[i][1] = 0.0;
  }
  fftw_execute(p_forward);

  G_hat = fftw_to_eigen(out);

  Eigen::VectorXcd F_hat(N);
  for (int i = 0; i < N; i++)
  {
    double k = i - N / 2;
    F_hat(i) = i_unit * k * G_hat(i);
    sum(i) = -F_hat(i) * exp(-i_unit * alpha * k * k * k * t);
  }
  for (int i = 0; i < N; i++)
  {
    int k_wave = i - N / 2;
    if (abs(k_wave) > N / 4)
      sum(i) = 0.0;
  }
  return sum;
}

int main()
{
  string filename, file2name;
  if (!is_linear_deleted)
  {
    filename = "K_dV_equation/data/fftw_data.csv";
    file2name = "K_dV_equation/data/fftw_data_time_step.csv";
  }
  else
  {
    filename = "K_dV_equation/data/fftw_linear_deleted_data.csv";
    file2name = "K_dV_equation/data/fftw_linear_deleted_data_time_step.csv";
  }
  ofstream file(filename);
  ofstream file2(file2name);

  in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);
  out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);
  p_forward = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  p_backward = fftw_plan_dft_1d(N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);

  Eigen::VectorXcd u(N);

  for (int i = 0; i < N; i++)
  {
    u(i) = sin(i * dx);
  }
  for (int i = 0; i < N; i++)
  {
    in[i][0] = u(i).real();
    in[i][1] = 0.0;
  }

  Eigen::VectorXcd u_hat(N);
  Eigen::VectorXcd u_record(N);
  Eigen::VectorXcd U_hat(N);
  fftw_execute(p_forward);
  u_hat = fftw_to_eigen(out);

  for (int i = 0; i < steps; i++)
  {
    double t = static_cast<double>(i) * dt;
    if (is_linear_deleted)
    {
      for (int j = 0; j < N; j++)
      {
        double k = j - N / 2;
        U_hat(j) = u_hat(j) * exp(-i_unit * alpha * k * k * k * t);
      }
      U_hat = runge_kutta_t_needed(t, U_hat, conversion_method_linear_delete);
      for (int l = 0; l < N; l++)
      {
        double k = l - N / 2;
        u_hat(l) = U_hat(l) * exp(i_unit * alpha * k * k * k * t);
      }
    }
    else
      u_hat = runge_kutta(u_hat, conversion_method);
    if (i % 1 == 0)
    {
      eigen_to_fftw(u_hat, in);
      fftw_execute(p_backward);

      for (int j = 0; j < N; j++)
      {
        u_record(j) = complex<double>(out[j][0], out[j][1]) / static_cast<double>(N);
        if (j != N - 1)
          file2 << u_record(j).real() << ",";
        else
          file2 << u_record(j).real();
      }
      file2 << endl;
    }
  }

  eigen_to_fftw(u_hat, in);
  fftw_execute(p_backward);

  for (int i = 0; i < N; i++)
  {
    u(i) = complex<double>(out[i][0], out[i][1]) / static_cast<double>(N);
  }

  for (int i = 0; i < N; i++)
  {
    file << i * dx << "," << u(i).real() << endl;
  }

  return 0;
}