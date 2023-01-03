
#include "1D_BTCS.h"

using namespace std;

void BTCS_1D()
{
	double x_l = -1.0;
	double x_r = 1.0;
	double dx = 0.025;
	int nx = ceil((x_r - x_l) / dx);
	dx = (x_r - x_l) / nx;

	double t = 1.0;
	double dt = 0.0025;
	int nt = ceil(t / dt);
	dt = t / nt;

	double alpha = 1 / (Pi * Pi);
	double CFL = alpha * dt / pow(dx, 2);
	vector<vector<double>> u(nt + 1, vector<double>(nx + 1, 0));// one timestep = one row
	vector<double> a(nx + 1, 0);
	vector<double> b(nx + 1, 0);
	vector<double> c(nx + 1, 0);
	vector<double> q(nx + 1, 0);

	// initial condition
	for (int i = 0; i < nx + 1; i++)
	{
		u[0][i] = -sin(Pi * (dx * i - 1));
	}
	a[0] = c[0] = q[0] = 0;
	b[0] = 1;
	for (int i = 1; i < nx; i++)
	{
		a[i] = c[i] = -CFL;
		b[i] = 1 + 2 * CFL;
	}
	a[nx] = c[nx] = q[nx] = 0;
	b[nx] = 1;
	for (int j = 1; j < nt + 1; j++)
	{
		for (int i = 1; i < nx; i++)
		{
			q[i] = u[j - 1][i];
		}
		thomasTridiagonal(a, b, c, u[j],q);
	}

	ofstream outfile("1d_BTCS_u.dat");
	if (outfile.is_open())
	{
		for (int i = 0; i < nx + 1; i++)
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