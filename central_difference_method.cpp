#include <iostream>
#include <numbers>
#include <vector>

using namespace std;

int N = 100;
double pi = std::numbers::pi;
int steps = 2000;
double dt = 2.0 * pi / steps;
double dx = 2.0 * pi / 2 / N;

// ∂u/∂t+∂u/∂x=0
int main()
{
  vector<double> u(2 * N + 1);
  vector<double> u_copy(2 * N + 1);
  vector<double> u_before(2 * N + 1);

  for (int i = 0; i < 2 * N + 1; i++)
  {
    u[i] = pow((1 - cos((i - N) * dx)) / 2, N);
  }

  for (int i = 0; i < steps; i++)
  {
    u[2 * N] = u[0];
    u_copy = u;
    for (int j = 0; j < 2 * N; j++)
    {
      int l = (j - 1 + 2 * N) % (2 * N);
      u[j] = u_copy[j] - dt * (u_copy[j + 1] - u_copy[l]) / 2 / dx;
    }
    u[2 * N] = u[0];
    u_before = u_copy;
    cout << u[0] << endl;
  }

  return 0;
}