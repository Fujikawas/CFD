#include <fftw3.h>
#include "1D_BTCS.h"
#include <complex>

void Poisson_FFT()
{
	double x_l = 0.0;
	double x_r = 1.0;
	int nx = 512;
	double dx = (x_r - x_l) / nx;
	vector<double> x(nx + 1, 0);
	for (int i = 0; i < nx + 1; i++)
	{
		x[i] = i* dx +x_l;
	}

	double y_b = 0.0;
	double y_t = 1.0;
	int ny = 512;
	double dy = (y_t - y_b) / ny;
	vector<double> y(ny + 1, 0);
	for (int i = 0; i < ny+1; i++)
	{
		y[i] = i * dy + y_b;
	}

	vector<vector<double>> ue(ny + 1, vector<double>(nx + 1, 0.0));
	vector<vector<double>> f(ny + 1, vector<double>(nx + 1, 0.0));
	vector<vector<double>> un(ny + 1, vector<double>(nx + 1, 0.0));
	double* fp = (double*)malloc((ny+1) * (nx+1) * sizeof(double));
	if (fp == nullptr) {
		std::cerr << "Memory allocation failed." << std::endl;
	}

	// analytic solution and initial condition
	double km = 16.0;
	double c1 = pow(1.0/km,2);
	double c2 = -8.0 * Pi * Pi;

	for (int i = 0; i < ny + 1; i++)
	{
		for (int j = 0; j < nx + 1; j++)
		{
			ue[i][j] = sin(2.0 * Pi * x[j]) * sin(2.0 * Pi * y[i]) + c1 * sin(km * 2.0 * Pi * x[j]) * sin(km * 2.0 * Pi * y[i]);
			f[i][j] = c2 * sin(2.0 * Pi * x[j]) * sin(2.0 * Pi * y[i]) + c2 * sin(km * 2.0 * Pi * x[j]) * sin(km * 2.0 * Pi * y[i]);
			fp[(nx+1) * i + j] = f[i][j];
		}
	}
	
	//FFT
	vector<double> kx(nx, 0);
	vector<double> ky(ny, 0);
	vector<vector<double>> e(ny, vector<double>(nx, 0.0));
	vector<vector<double>> u(ny, vector<double>(nx, 0.0));
	double aa = -2.0 / dx / dx - 2.0 / dy / dy;
	double bb = 2.0 / dx / dx;
	double cc = 2.0 / dy / dy;
	double hx = 2.0 * Pi / nx;
	double eps = 1.0e-6;

	
	for (int i = 0; i < nx/2; i++)
	{
		kx[i] = hx * i;
		kx[i + nx / 2] = hx * (i - nx / 2);
	}
	kx[0] = eps;
	ky = kx;

	//FFT
	fftw_complex* in, * out;
	fftw_plan p;
	in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * nx * ny);
	
	for (int i = 0; i < ny; i++)
	{
		for (int j = 0; j < nx; j++)
		{
			in[nx * i + j][0] = fp[(nx + 1) * i + j];
			in[nx * i + j][1] = 0.0;
		}
	}
	free(fp);
	out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * nx * ny);
	p = fftw_plan_dft_2d(ny, nx, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(p); 
	out[0][0] = 0.0;
	for (int i = 0; i < ny; i++)
	{
		for (int j = 0; j < nx; j++)
		{
			out[nx * i + j][0] /= aa + bb * cos(kx[j]) + cc * cos(ky[i]);
			out[nx * i + j][1] /= aa + bb * cos(kx[j]) + cc * cos(ky[i]);
			
		}
	}
	fftw_plan pi;
	pi = fftw_plan_dft_2d(ny, nx, out, in, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(pi);
	fftw_destroy_plan(p);
	fftw_destroy_plan(pi);
	for (int i = 0; i < ny; i++)
	{
		for (int j = 0; j < nx; j++)
		{
			un[i][j]= in[nx * i + j][0]/nx/ny;
		}
		un[i][nx] = un[i][0];
	}
	for (int i = 0; i < nx+1; i++)
	{
		un[ny][i] = un[0][i];
	}

	//write
	ofstream outfile("Poisson_FFT.dat");
	if (outfile.is_open())
	{
		for (int i = 0; i < ny + 1; i++)
		{
			for (int j = 0; j < nx + 1; j++)
			{
				outfile << un[j][i] << " ";
			}
			outfile << endl;
		}
		outfile << endl;
	}
	else
	{
		std::cerr << "Error: unable to open file for writing" << std::endl;
	}

	fftw_free(in); fftw_free(out);

	return;
}
