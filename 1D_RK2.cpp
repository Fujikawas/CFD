#include "1D_BTCS.h"
using namespace std;

void RK2_1D()
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
	vector<vector<double>> u_half(nt + 1, vector<double>(nx + 1, 0));// one timestep = one row

	// initial condition
	for (int i = 0; i < nx + 1; i++)
	{
		u[0][i] = -sin(Pi * (dx * i - 1));
	}
	for (int j = 1; j < nt + 1; j++)// one time step
	{
		u_half[j][0] = 0.0; //BC
		for (int i = 1; i < nx; i++)
		{
			u_half[j][i] = u[j - 1][i] + CFL / 2 * (u[j - 1][i + 1] - 2 * u[j - 1][i] + u[j - 1][i - 1]);
		}
		u_half[j][nx] = 0.0; //BC
		u[j][0] = 0.0;
		for (int i = 1; i < nx; i++)
		{
			u[j][i] = u[j - 1][i] + CFL * (u_half[j][i + 1] - 2 * u_half[j][i] + u_half[j][i - 1]);
		}
		u[j][nx] = 0.0;

	}

	ofstream outfile("1d_RK2_u.dat");
	if (outfile.is_open())
	{
		/*vector<vector<double>>::iterator ptr;
		for (ptr = u.begin(); ptr < u.end(); ptr++)
		{
			ostream_iterator<double> output_iterator(outfile, "\n");
			copy((*ptr).begin(), (*ptr).end(), output_iterator);
		}
		outfile.close();*/
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