#include <Eigen/Dense>
#include <fftw3.h>
#include "constant.h"
#include <complex>

fftw_complex *in, *out;
fftw_plan p_forward, p_backward;

Eigen::VectorXcd fftw_to_eigen(fftw_complex *out)
{
  Eigen::VectorXcd result(M);
  for (int i = 0; i < M; i++)
  {
    int k_wave = i - M / 2;
    int fft_idx = (k_wave >= 0) ? k_wave : (k_wave + M);
    result(i) = complex<double>(out[fft_idx][0], out[fft_idx][1]);
  }
  return result;
}

void eigen_to_fftw(Eigen::VectorXcd &u_hat, fftw_complex *input)
{
  for (int i = 0; i < M; i++)
  {
    int k_wave = i - M / 2;
    int fft_idx = (k_wave >= 0) ? k_wave : (k_wave + M);
    input[fft_idx][0] = u_hat(i).real();
    input[fft_idx][1] = u_hat(i).imag();
  }
}