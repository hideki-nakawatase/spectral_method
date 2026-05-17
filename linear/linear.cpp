#include <iostream>
#include <numbers>
#include <vector>
#include <fstream>

using namespace std;

int N = 100;
double pi = std::numbers::pi;
int steps = 200000;
double dt = 2.0 * pi / steps;
double dx = 2.0 * pi / (2 * N + 1);

// ∂u/∂t+∂u/∂x=0
int main()
{
  ofstream ofs("linear/linear.csv");
  vector<double> u(2 * N + 1), u_true(2 * N + 1), u_copy(2 * N + 1);

  for (int i = 0; i < 2 * N + 1; i++)
  {
    u[i] = pow((1 - cos(i * dx)) / 2, N);
    u_true[i] = u[i];
  }

  for (int i = 0; i < steps; i++)
  {
    u_copy = u;
    for (int j = 0; j < 2 * N + 1; j++)
    {
      int l = (j - 1 + 2 * N + 1) % (2 * N + 1);
      int k = (j + 1) % (2 * N + 1);
      u[j] -= dt * (u_copy[k] - u_copy[l]) / 2 / dx;
    }
  }
  for (int i = 0; i < 2 * N + 1; i++)
  {
    ofs << i * dx << "," << u[i] << "," << u_true[i] << endl;
  }
  return 0;
}