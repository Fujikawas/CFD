#include <fftw3.h>
#include "1D_BTCS.h"

void Poisson_MG()
{
	double x_l = 0.0;
	double x_r = 1.0;
	int nx = 128;
	double dx = (x_r - x_l) / nx;
	vector<double> x(nx + 1, 0);
	for (int i = 0; i < nx + 1; i++)
	{
		x[i] = i * dx + x_l;
	}

	double y_b = 0.0;
	double y_t = 1.0;
	int ny = 128;
	double dy = (y_t - y_b) / ny;
	vector<double> y(ny + 1, 0);
	for (int i = 0; i < ny + 1; i++)
	{
		y[i] = i * dy + y_b;
	}

	double tolerance = 1.0e-4;
	int max_iter = 10000;
	int v1 = 2;//relaxation
	int v2 = 2;//prolongation
	int v3 = 2;//coarsest level

	vector<vector<double>> ue(ny + 1, vector<double>(nx + 1, 0.0));
	vector<vector<double>> f(ny + 1, vector<double>(nx + 1, 0.0));
	vector<vector<double>> un(ny + 1, vector<double>(nx + 1, 0.0));


	// analytic solution and initial condition
	double km = 16.0;
	double c1 = pow(1.0 / km, 2);
	double c2 = -2.0 * Pi * Pi;

	for (int i = 0; i < ny + 1; i++)
	{
		for (int j = 0; j < nx + 1; j++)
		{
			ue[i][j] = sin(2.0 * Pi * x[j]) * sin(2.0 * Pi * y[i]) + c1 * sin(km * Pi * x[j]) * sin(km * Pi * y[i]);
			f[i][j] = 4.0 * c2 * sin(2.0 * Pi * x[j]) * sin(2.0 * Pi * y[i]) + c2 * sin(km * Pi * x[j]) * sin(km * Pi * y[i]);

		}
	}
	for (int i = 0; i < ny + 1; i++)
	{
		un[i][0] = ue[i][0];
		un[i][nx] = ue[i][nx];
	}
	for (int i = 0; i < nx + 1; i++)
	{
		un[0][i] = ue[0][i];
		un[ny][i] = ue[ny][i];
	}

	vector<vector<double>> r(ny + 1, vector<double>(nx + 1, 0.0));
	double init_rms = 0.0;
	double rms = 0.0;

	for (int i = 1; i < ny; i++)
	{
		for (int j = 1; j < nx; j++)
		{
			double d2udy2 = (un[i + 1][j] - 2 * un[i][j] + un[i - 1][j]) / dy / dy;
			double d2udx2 = (un[i][j + 1] - 2 * un[i][j] + un[i][j - 1]) / dx / dx;
			r[i][j] = f[i][j] - d2udx2 - d2udy2;
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

	//allocate memory for grid size at different levels
	vector<int> lnx(2, 0); 
	vector<int> lny(2, 0); 
	vector<double> ldx(2, 0.0); 
	vector<double> ldy(2, 0.0); 

	lnx[0] = nx;
	lny[0] = ny;
	lnx[1] = lnx[0] / 2;
	lny[1] = lny[0] / 2;
	ldx[0] = dx;
	ldy[0] = dy;
	ldx[1] = ldx[0] * 2;
	ldy[1] = ldy[0] * 2;

	//allocate matrix for storage at fine level, residual at fine level is already defined at global level
	vector<vector<double>> prol_fine(lny[0] + 1, vector<double>(lnx[0] + 1, 0.0));

	//allocate matrix for storage at coarse levels
	vector<vector<double>> fc(lny[1] + 1, vector<double>(lnx[1] + 1, 0.0));
	vector<vector<double>> unc(lny[1] + 1, vector<double>(lnx[1] + 1, 0.0));


	int iter_count = 0;
	double exp_rms = tolerance * init_rms;

	for (iter_count = 0; iter_count < max_iter && rms>exp_rms; iter_count++)
	{
		//call relaxation on fine grid and compute the numerical solution for fixed number of iteration
		GaussSeidel_MG(lnx[0], lny[0], dx, dy, f, un, v1);
		//compute residual
		rms = 0.0;
		for (int i = 1; i < ny; i++)
		{
			for (int j = 1; j < nx; j++)
			{
				double d2udy2 = (un[i + 1][j] - 2 * un[i][j] + un[i - 1][j]) / dy / dy;
				double d2udx2 = (un[i][j + 1] - 2 * un[i][j] + un[i][j - 1]) / dx / dx;
				r[i][j] = f[i][j] - d2udx2 - d2udy2;

			}
		}
		for (int i = 1; i < ny; i++)
		{
			for (int j = 1; j < nx; j++)
			{
				rms += r[i][j] * r[i][j];
			}
		}
		rms = sqrt(rms / (nx - 1) / (ny - 1));

		//restrict the residual from fine level to coarse level
		restriction(lnx[0], lny[0], lnx[1], lny[1], r, fc);
		//solve on the coarsest level and relax V3 times
		GaussSeidel_MG(lnx[1], lny[1], ldx[1], ldy[1], fc, unc, v3);
		//prolongate solution from coarse level to fine level
		prolongation(lnx[1], lny[1], lnx[0], lny[0], unc, prol_fine);
		//correct the solution on fine level
		for (int i = 1; i < lny[0]; i++)
		{
			for (int j = 1; j < lnx[0]; j++)
			{
				un[i][j] += prol_fine[i][j];
			}
		}
		//relax v2 times
		GaussSeidel_MG(lnx[0], lny[0], dx, dy, f, un, v2);

		cout << "iteration times " << iter_count << endl;
		cout << "residual " << rms << endl;
	}
	cout << "iteration times until convergence:" << iter_count << endl;

	//write
	ofstream outfile("Poisson_MG.dat");
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

void restriction(int nxf, int nyf, int nxc, int nyc, vector<vector<double>> r, vector<vector<double>>& ec)
{
	for (int i = 1; i < nyc; i++)
	{
		for (int j = 1; j < nxc; j++)
		{
			double center = 4.0 * r[2 * i][2 * j];
			double grid = 2.0 * (r[2 * i][2 * j + 1] + r[2 * i][2 * j - 1] + r[2 * i + 1][2 * j] + r[2 * i - 1][2 * j]);
			double corner = 1.0 * (r[2 * i +1][2 * j + 1] + r[2 * i +1][2 * j - 1] + r[2 * i - 1][2 * j + 1] + r[2 * i - 1][2 * j - 1]);
			ec[i][j] = (center + grid + corner) / 16.0;
		}
	}
	for (int i = 0; i < nxc + 1; i++)
	{
		ec[0][i] = r[0][2 * i];
		ec[nyc][i] = r[nyf][2 * i];
	}
	for (int i = 0; i < nyc + 1; i++)
	{
		ec[i][0] = r[2 * i][0];
		ec[i][nxc] = r[2 * i][nxf];
	}
	return;
}
void prolongation(int nxc, int nyc, int nxf, int nyf, vector<vector<double>> unc, vector<vector<double>> ef)
{
	for (int i = 0; i < nyc; i++)
	{
		for (int j = 0; j < nxc; j++)
		{
			ef[2 * i][2 * j] = unc[i][j];
			ef[2 * i][2 * j + 1] = 0.5 * (unc[i][j] + unc[i][j + 1]);
			ef[2 * i + 1][2 * j] = 0.5 * (unc[i][j] + unc[i + 1][j]);
			ef[2 * i + 1][2 * j + 1] = 0.25 * (unc[i][j] + unc[i][j + 1] + unc[i + 1][j] + unc[i + 1][j + 1]);
		}
	}

	for (int i = 0; i < nxc + 1; i++)
	{
		ef[2 * i][0] = unc[i][0];
		ef[2 * i][nyf] = unc[i][nxc];
	}
	for (int i = 0; i < nxc + 1; i++)
	{
		ef[0][2 * i] = unc[0][i];
		ef[nyf][2 * i] = unc[nyc][i];
	}
	return;
}
void GaussSeidel_MG(int nx, int ny, double dx, double dy, vector<vector<double>> f, vector<vector<double>>& un, int V)
{
	vector<vector<double>> rt(ny + 1, vector<double>(nx + 1, 0.0));
	double den = -2.0 / dx / dx - 2.0 / dy / dy;
	for (int iter = 0; iter < V; iter++)
	{
		//compute solution at next time step 
		for (int i = 1; i < ny; i++)
		{
			for (int j = 1; j < nx; j++)
			{
				rt[i][j] = f[i][j] - (un[i + 1][j] - 2 * un[i][j] + un[i - 1][j]) / dy / dy - (un[i][j + 1] - 2 * un[i][j] + un[i][j - 1]) / dx / dx;
				un[i][j] += rt[i][j] / den;
			}
		}
	}
	return;
}