#include "constant.h"
#include "fourier.h"
#include "runge_kutta.h"
#include <fstream>
#include <Eigen/Dense>
#include <iostream>
#include <chrono>
#include <fftw3.h>

using namespace std;

bool is_fft = true;

fftw_complex *in, *out;
fftw_plan p_forward, p_backward;

Eigen::VectorXcd fftw_to_eigen(fftw_complex *out)
{
    Eigen::VectorXcd result(M);
    for (int i = 0; i < M; i++)
    {
        int k_wave = i - N;
        int fft_idx = (k_wave >= 0) ? k_wave : (k_wave + M);
        result(i) = complex<double>(out[fft_idx][0], out[fft_idx][1]);
    }
    return result;
}

void eigen_to_fftw(Eigen::VectorXcd &u_hat, fftw_complex *input)
{
    for (int i = 0; i < M; i++)
    {
        int k_wave = i - N;
        int fft_idx = (k_wave >= 0) ? k_wave : (k_wave + M);
        input[fft_idx][0] = u_hat(i).real();
        input[fft_idx][1] = u_hat(i).imag();
    }
}

Eigen::VectorXcd conversion_method(Eigen::VectorXcd &u_hat)
{
    Eigen::VectorXcd G(M), G_hat(M), u(M);
    if (is_fft)
    {
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
    }
    else
    {
        u = inverse_fourier(u_hat);
        for (int i = 0; i < M; i++)
        {
            G(i) = (u(i).real() * u(i).real()) / 2.0;
        }
        G_hat = forward_fourier(G);
    }
    Eigen::VectorXcd F_hat(M);
    for (int i = 0; i < M; i++)
    {
        int k_wave = i - N;
        F_hat(i) = i_unit * static_cast<double>(k_wave) * G_hat(i);
    }
    return -F_hat;
}

int main()
{
    in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * M);
    out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * M);

    p_forward = fftw_plan_dft_1d(M, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    p_backward = fftw_plan_dft_1d(M, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    auto start = chrono::steady_clock::now();

    ofstream file(is_fft ? "non_linear/conversion_method_fftw.csv" : "non_linear/conversion_method.csv");

    if (!file.is_open())
        cout << "file not found" << endl;

    Eigen::VectorXcd u(M);
    u.setZero();
    for (int i = 0; i < M; i++)
        u(i) = sin(i * dx);

    Eigen::VectorXcd u_hat(M);
    if (is_fft)
    {
        for (int i = 0; i < M; i++)
        {
            in[i][0] = u(i).real();
            in[i][1] = u(i).imag();
        }
        fftw_execute(p_forward);
        u_hat = fftw_to_eigen(out);
    }
    else
    {
        u_hat = forward_fourier(u);
    }

    for (int i = 0; i < steps; i++)
    {
        u_hat = runge_kutta(u_hat, conversion_method);
    }

    if (is_fft)
    {
        eigen_to_fftw(u_hat, in);
        fftw_execute(p_backward);
        fftw_to_eigen(out);
        for (int i = 0; i < M; i++)
        {
            u(i) = complex<double>(out[i][0], out[i][1]) / static_cast<double>(M);
        }
    }
    else
    {
        u = inverse_fourier(u_hat);
    }

    auto end = chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    cout << duration / 1000000 << endl;

    for (int i = 0; i < M; i++)
    {
        file << i * dx << "," << u(i).real() << endl;
    }

    return 0;
}