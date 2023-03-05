#include <fftw3.h>
#include "1D_BTCS.h"

void Poisson_CG()
{
	double x_l = 0.0;
	double x_r = 1.0;
	int nx = 512;
	double dx = (x_r - x_l) / nx;
	vector<double> x(nx + 1, 0);
	for (int i = 0; i < nx + 1; i++)
	{
		x[i] = i * dx + x_l;
	}

	double y_b = 0.0;
	double y_t = 1.0;
	int ny = 512;
	double dy = (y_t - y_b) / ny;
	vector<double> y(ny + 1, 0);
	for (int i = 0; i < ny + 1; i++)
	{
		y[i] = i * dy + y_b;
	}

	double tolerance = 1.0e-4;
	int max_iter = 10000;
	double tiny = 1.0e-16;

	vector<vector<double>> ue(ny + 1, vector<double>(nx + 1, 0.0));
	vector<vector<double>> f(ny + 1, vector<double>(nx + 1, 0.0));
	vector<vector<double>> un(ny + 1, vector<double>(nx + 1, 0.0));


	// analytic solution and initial condition
	for (int i = 0; i < ny + 1; i++)
	{
		for (int j = 0; j < nx + 1; j++)
		{
			ue[i][j] = (x[j] * x[j] - 1.0) * (y[i] * y[i] - 1.0);
			f[i][j] = -2.0 * (2.0 - x[j] * x[j] - y[i] * y[i]);
		}
	}
	for (int i = 0; i < ny + 1; i++)
	{
		un[i][0] = ue[i][0];
		un[i][ny] = ue[i][ny];
	}
	for (int i = 0; i < nx + 1; i++)
	{
		un[0][i] = ue[0][i];
		un[nx][i] = ue[nx][i];
	}

	vector<vector<double>> r(ny + 1, vector<double>(nx + 1, 0.0));
	vector<vector<double>> p(ny + 1, vector<double>(nx + 1, 0.0));//matrix for direction
	vector<vector<double>> del_p(ny + 1, vector<double>(nx + 1, 0.0));//matrix for direction
	double init_rms = 0.0;
	double rms = 0.0;

	for (int i = 1; i < ny; i++)
	{
		for (int j = 1; j < nx; j++)
		{
			double d2udx2 = (un[i + 1][j] - 2 * un[i][j] + un[i - 1][j]) / dx / dx;
			double d2udy2 = (un[i][j + 1] - 2 * un[i][j] + un[i][j - 1]) / dy / dy;
			r[i][j] = f[i][j] - d2udx2 - d2udy2;
			p[i][j] = r[i][j];

		}
	}
	//compute residual
	for (int i = 1; i < ny; i++)
	{
		for (int j = 1; j < nx; j++)
		{
			init_rms += r[i][j] * r[i][j];
		}
	}
	init_rms = sqrt(init_rms / (nx - 1) / (ny - 1));
	rms = init_rms;

	int iter_count = 0;
	double den = -2.0 / dx / dx - 2.0 / dy / dy;
	double exp_rms = tolerance * init_rms;
	
	for (iter_count = 0; iter_count < max_iter && rms>exp_rms; iter_count++)
	{
		//correct solution
		for (int i = 1; i < ny; i++)
		{
			for (int j = 1; j < nx; j++)
			{
				del_p[i][j] = (p[i + 1][j] - 2 * p[i][j] + p[i - 1][j]) / dx / dx + (p[i][j + 1] - 2 * p[i][j] + p[i][j - 1]) / dy / dy;
			}
		}
		//compute new residual
		double aa = 0.0;
		double bb = 0.0;
		for (int i = 1; i < ny; i++)
		{
			for (int j = 1; j < nx; j++)
			{
				aa += r[i][j] * r[i][j];
				bb += del_p[i][j] * p[i][j];
			}
		}
		double cc = aa / (bb + tiny);
		//update the numerical solution by adding some component of conjugate vector
		for (int i = 1; i < ny; i++)
		{
			for (int j = 1; j < nx; j++)
			{
				un[i][j] += cc * p[i][j];
			}
		}
		bb = aa;
		aa = 0.0;
		//update the residual by removing some component of previous residual
		for (int i = 1; i < ny; i++)
		{
			for (int j = 1; j < nx; j++)
			{
				r[i][j] -= cc * del_p[i][j];
				aa += r[i][j] * r[i][j];
			}
		}
		cc = aa / (bb + tiny);
		//update the conjugate vector
		for (int i = 0; i < ny; i++)
		{
			for (int j = 0; j < nx; j++)
			{
				p[i][j] = r[i][j] + cc * p[i][j];
			}
		}
		rms = 0.0;
		for (int i = 1; i < ny; i++)
		{
			for (int j = 1; j < nx; j++)
			{

				rms += r[i][j] * r[i][j];
			}
		}
		rms = sqrt(rms / (nx - 1) / (ny - 1));
		cout << "iteration times " << iter_count << endl;
		cout << "residual " << rms << endl;
	}
	cout << "iteration times until convergence:" << iter_count << endl;

	//write
	ofstream outfile("Poisson_CG.dat");
	if (outfile.is_open())
	{
		for (int i = 0; i < ny + 1; i++)
		{
			for (int j = 0; j < nx + 1; j++)
			{
				outfile << un[i][j] << " ";
			}
			outfile << endl;
		}
		outfile << endl;
	}
	else
	{
		std::cerr << "Error: unable to open file for writing" << std::endl;
	}
	return;
}
