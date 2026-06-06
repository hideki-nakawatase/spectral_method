#include <numbers>
#include <Eigen/Dense>
#include <fftw3.h>
#include <complex>
#include <conversion_fftw_eigen.h>

using namespace std;

fftw_complex *in, *out;
fftw_plan p_forward, p_backward;

double L = 2 * numbers::pi;
int N = 256;
int k_wave_lim = 85;
int steps = 4000;
int dx = L / N;
double alpha = 0.01;
complex<double> i_unit(0.0, 1.0);

Eigen::VectorXcd conversion_method(Eigen::VectorXcd &u_hat)
{
  Eigen::VectorXcd G(M), G_hat(M), u(M), sum(M);
  sum.setZero();
  eigen_to_fftw(u_hat, in);
  fftw_execute(p_backward);

  for (int i = 0; i < M; i++)
  {
    u(i) = complex<double>(out[i][0], out[i][1]) / static_cast<double>(M);
  }
  for (int i = 0; i < M; i++)
  {
    in[i][0] = u(i).real() * u(i).real() / 2.0;
    in[i][1] = 0.0;
  }
  fftw_execute(p_forward);

  G_hat = fftw_to_eigen(out);

  Eigen::VectorXcd F_hat(M);
  for (int i = 0; i < M; i++)
  {
    int k_wave = i - N;
    double k = static_cast<double>(k_wave);
    F_hat(i) = i_unit * k * G_hat(i);
    sum(i) += F_hat(i) - alpha * i_unit * k * k * k * u_hat(i);
  }
  return -sum;
}

int main()
{
  in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);
  out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N);
  p_forward = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  p_backward = fftw_plan_dft_1d(N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
  Eigen::VectorXcd u(N);
  u.setZero();
  for (int i = 0; i < N; i++)
  {
    u(i) = sin(i * dx);
  }

  Eigen::VectorXcd u_hat(N);
  for (int i = 0; i < N; i++)
  {
    in[i][0] = u(i).real();
    in[i][1] = u(i).real();
  }

  fftw_execute(p_forward);
  u_hat = fftw_to_eigen(out);

  return 0;
}