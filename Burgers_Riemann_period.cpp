#include "1D_BTCS.h"

using namespace std;

void Riemann_period()
{
	double x_l = 0;
	double x_r = 1.0;
	int nx = 200;
	double dx = (x_r - x_l) / nx;
	vector<double> x(nx, 0);//temperory array by RK3 scheme

	double t = 0.25;
	double dt = 0.0001;
	int nt = ceil(t / dt);
	dt = t / nt;
	vector<vector<double>> u(nt + 1, vector<double>(nx, 0));// one timestep = one row, FDM has nx+1 columns
	vector<double> ut(nx, 0);//temperory array by RK3 scheme
	vector<double> r(nx, 0);//general spatial FD 

	for (int i = 0; i < nx; i++)
	{
		x[i] = (i + 0.5) * dx;
		u[0][i] = sin(2.0 * Pi * x[i]);
	}

	for (int j = 1; j < nt + 1; j++)// one time step
	{
		rhs_RM(nx, dx, u[j - 1], r);
		for (int i = 1; i < nx - 1; i++)
		{
			ut[i] = u[j - 1][i] + dt * r[i];
		}
		rhs_RM(nx, dx, ut, r);
		for (int i = 1; i < nx - 1; i++)
		{
			ut[i] = 0.75 * u[j - 1][i] + 0.25 * ut[i] + dt / 4.0 * r[i];
		}
		rhs_RM(nx, dx, ut, r);
		for (int i = 1; i < nx - 1; i++)
		{
			u[j][i] = u[j - 1][i] / 3.0 + 2.0 / 3.0 * ut[i] + dt * 2.0 / 3.0 * r[i];
		}

	}

	ofstream outfile("Riemann_period.dat");
	if (outfile.is_open())
	{
		for (int i = 0; i < nx; i++)
		{
			outfile << u[nt][i] << " ";
		}
		outfile << endl;
	}
	else
	{
		std::cerr << "Error: unable to open file for writing" << std::endl;
	}
	return;
}

void rhs_RM(int nx, double dx, vector<double> u, vector<double>& r)
{
	vector<double> uL(nx + 1, 0);
	vector<double> uR(nx + 1, 0);
	//left and right side fluxes at the interface
	vector<double> fL(nx + 1, 0);
	vector<double> fR(nx + 1, 0);
	vector<double> f(nx + 1, 0);

	uL = wenoL_FS(nx, u);
	uR = wenoR_FS(nx, u);

	for (int i = 0; i < nx + 1; i++)
	{
		fL[i] = 0.5 * uL[i] * uL[i];
		fR[i] = 0.5 * uR[i] * uR[i];
	}
	
	Riemann(nx,u,uL,uR,f,fL,fR);

	for (int i = 0; i < nx; i++)
	{
		r[i] = -(f[i + 1] - f[i]) / dx;
	}
	return;
}

void Riemann(int nx, vector<double> u, vector<double> uL, vector<double> uR, vector<double> &f, vector<double> fL, vector<double> fR)
{
	vector<double> c(nx + 1, 0);
	for (int i = 1; i < nx; i++)
	{
		c[i] = max(abs(u[i - 1]), abs(u[i]));
	}
	c[0] = max(abs(u[0]), abs(u[nx - 1]));
	c[nx] = max(abs(u[0]), abs(u[nx - 1]));
	
	//Interface fluxes(Rusanov)
	for (int i = 0; i < nx + 1; i++)
	{
		f[i] = 0.5 * (fR[i] + fL[i]) - 0.5 * c[i] * (uR[i] - uL[i]);
	}
}