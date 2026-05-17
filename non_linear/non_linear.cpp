#include <iostream>
#include <numbers>
#include <vector>
#include <fstream>
#include <algorithm>
#include <complex>

using namespace std;

int N = 100;
int steps = 200000;
double dt = 1.0 / steps;
double dx = 2.0 / (2 * N + 1);
std::complex<double> i_unit = 1.0i;

vector<complex<double>> wave(int k)
{
  vector<complex<double>> e(2 * N + 1);
  for (int i = 0; i < 2 * N + 1; i++)
    e[i] = cos(k * i * dx) + i_unit * sin(k * i * dx);
  return e;
}

// ∂u/∂t+u*∂u/∂x=0

void spectral_method(ofstream ofs) {}

int main()
{
  ofstream ofs("non_linear/non_linear.csv");
  vector<complex<double>> u(2 * N + 1);
  vector<complex<double>> u_hat(2 * N + 1), u_hat_copy(2 * N + 1), e(2 * N + 1);

  for (int i = 0; i < 2 * N + 1; i++)
  {
    u[i] = sin(i + dx);
  }

  for (int i = 0; i < steps; i++)
  {
    std::complex<double> sum;
    for (int j = 0; j < 2 * N + 1; j++)
    {
      u_hat_copy = u_hat;
      int low_lm = max(0, j - N);
      int up_lm = min(2 * N, j + N);
      std::complex<double> sum = 0;
      for (int k = low_lm; k < up_lm; k++)
      {
        sum += i_unit * static_cast<double>(k) * u_hat_copy[j - k] * u_hat_copy[k];
      }
      u_hat[j] -= dt * sum;
    }
  }

  for (int i = 0; i < 2 * N + 1; i++)
    u[i] = 0.0;

  for (int j = 0; j < 2 * N + 1; j++)
  {
    vector<complex<double>> e_wave = wave(j);
    for (int i = 0; i < 2 * N + 1; i++)
      u[i] += u_hat[j] * e_wave[i];
  }

  for (int i = 0; i < 2 * N + 1; i++)
  {
    ofs << i * dx << "," << u[i] << endl;
  }
  return 0;
}