#include <Eigen/Dense>
#include <fstream>
#include "constant.h"
#include "runge_kutta.h"

using namespace std;

Eigen::VectorXd center_diff(Eigen::VectorXd &u)
{
    Eigen::VectorXd du(M);

    for (int i = 0; i < M; i++)
    {
        int p1 = (i + 1) % M;
        int m1 = (i - 1 + M) % M;
        du[i] = -u(i) * (u(p1) - u(m1)) / 2 / dx;
    }
    return du;
}

int main()
{
    ofstream file("non_linear/center_diff.csv");

    Eigen::VectorXd u(M), du1(M), du2(M), du3(M), du4(M), k1(M), k2(M), k3(M), k4(M), tmp(M);

    for (int i = 0; i < M; i++)
        u(i) = sin(i * dx);

    for (int i = 0; i < steps; i++)
    {
        u = runge_kutta(u, center_diff);
    }

    for (int i = 0; i < M; i++)
        file << i * dx << "," << u(i) << "\n";

    return 0;
}