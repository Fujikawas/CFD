#include "1D_BTCS.h"

using namespace std;

void WENO5_Dirichlet()
{
	double x_l = 0;
	double x_r = 1.0;
	int nx = 200;
	double dx = (x_r - x_l) / nx;

	double t = 0.25;
	double dt = 0.0001;
	int nt = ceil(t / dt);
	dt = t / nt;

	vector<vector<double>> u(nt + 1, vector<double>(nx + 1, 0));// one timestep = one row
	vector<double> ut(nx + 1, 0);//temperory array by RK3 scheme
	vector<double> r(nx + 1, 0);//general spatial FD 

	// initial condition
	for (int i = 0; i < nx + 1; i++)
	{
		u[0][i] = sin(2 * Pi * dx * i);
	}

	for (int j = 1; j < nt + 1; j++)// one time step
	{
		u[j][0] = 0.0; //BC
		u[j][nx] = 0.0; //BC
		
		rhs(nx, dx, u[j - 1], r);		
		for (int i = 1; i < nx; i++)
		{
			ut[i] = u[j - 1][i] + dt * r[i];
		}
		rhs(nx, dx, ut, r);
		for (int i = 1; i < nx; i++)
		{
			ut[i] = 0.75 * u[j - 1][i] + 0.25 * ut[i] + dt / 4.0 * r[i];
		}
		rhs(nx, dx, ut, r);		
		for (int i = 1; i < nx; i++)
		{
			u[j][i] = u[j - 1][i] / 3.0 + 2.0 / 3.0 * ut[i] + dt * 2.0 / 3.0 * r[i];
		}

	}

	ofstream outfile("WENO5.dat");
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

void rhs(int nx, double dx, vector<double> u, vector<double> &r)
{
	vector<double> uL(nx, 0);
	vector<double> uR(nx+1, 0);

	wenoL(nx, u, uL);
	wenoR(nx, u, uR);

	for (int i = 1; i < nx; i++)
	{
		if (u[i] >= 0.0)
			r[i] = -u[i] * (uL[i] - uL[i - 1]) / dx;
		else
			r[i] = -u[i] * (uR[i + 1] - uR[i]) / dx;
	}
	return;
}

void wenoL(int nx, vector<double> u, vector<double> &uL)
{
	int i = 0;
	double	v1 = 3.0 * u[i] - 2.0 * u[i + 1];
	double 	v2 = 2.0 * u[i] - u[i + 1];
	double	v3 = u[i];
	double	v4 = u[i + 1];
	double 	v5 = u[i + 2];
	uL[i] = wcL(v1, v2, v3, v4, v5);

	i = 1;
	v1 = 2.0 * u[i - 1] - u[i];
	v2 = u[i - 1];
	v3 = u[i];
	v4 = u[i + 1];
	v5 = u[i + 2];
	uL[i] = wcL(v1, v2, v3, v4, v5);

	for (i = 2; i < nx-1;i++)
	{
		v1 = u[i - 2];
		v2 = u[i - 1];
		v3 = u[i];
		v4 = u[i + 1];
		v5 = u[i + 2];
		uL[i] = wcL(v1, v2, v3, v4, v5);
	}
	i = nx - 1;
	v1 = u[i - 2];
	v2 = u[i - 1];
	v3 = u[i];
	v4 = u[i + 1];
	v5 = 2.0 * u[i + 1] - u[i];
	uL[i] = wcL(v1, v2, v3, v4, v5);

	return;
}
void wenoR(int nx, vector<double> u, vector<double> &uR)
{
	int i = 1;
	double v1 = 2.0 * u[i - 1] - u[i];
	double v2 = u[i - 1];
	double v3 = u[i];
	double v4 = u[i + 1];
	double v5 = u[i + 2];
	uR[i] = wcR(v1, v2, v3, v4, v5);

	for (i = 2; i < nx-1; i++)
	{
		v1 = u[i - 2];
		v2 = u[i - 1];
		v3 = u[i];
		v4 = u[i + 1];
		v5 = u[i + 2];
		uR[i] = wcR(v1, v2, v3, v4, v5);
	}
	i = nx-1;
	v1 = u[i - 2];
	v2 = u[i - 1];
	v3 = u[i];
	v4 = u[i + 1];
	v5 = 2.0 * u[i + 1] - u[i];
	uR[i] = wcR(v1, v2, v3, v4, v5);

	i = nx;
	v1 = u[i - 2];
	v2 = u[i - 1];
	v3 = u[i];
	v4 = 2.0 * u[i] - u[i - 1];
	v5 = 3.0 * u[i] - 2.0 * u[i - 1];
	uR[i] = wcR(v1, v2, v3, v4, v5);

	return;
}

double wcL(double v1, double v2, double v3, double v4, double v5)
{
	double eps = 1.0e-6;

		//smoothness indicators
	double s1 = 13.0 / 12.0 * pow(v1 - 2.0 * v2 + v3,2) + 0.25 * pow(v1 - 4.0 * v2 + 3.0 * v3,2);
	double s2 = 13.0 / 12.0 * pow(v2 - 2.0 * v3 + v4,2) + 0.25 * pow(v2 - v4,2);
	double s3 = 13.0 / 12.0 * pow(v3 - 2.0 * v4 + v5,2) + 0.25 * pow(3.0 * v3 - 4.0 * v4 + v5,2);

		//computing nonlinear weights w1, w2, w3
	double c1 = 1.0e-1 / (pow(eps + s1,2));
	double c2 = 6.0e-1 / (pow(eps + s2,2));
	double c3 = 3.0e-1 / (pow(eps + s3,2));

	double w1 = c1 / (c1 + c2 + c3);
	double w2 = c2 / (c1 + c2 + c3);
	double w3 = c3 / (c1 + c2 + c3);

		//candiate stencils
	double q1 = v1 / 3.0 - 7.0 / 6.0 * v2 + 11.0 / 6.0 * v3;
	double q2 = -v2 / 6.0 + 5.0 / 6.0 * v3 + v4 / 3.0;
	double q3 = v3 / 3.0 + 5.0 / 6.0 * v4 - v5 / 6.0;

		//reconstructed value at interface
	double f = w1 * q1 + w2 * q2 + w3 * q3;

	return f;
}
double wcR(double v1, double v2, double v3, double v4, double v5)
{
	double eps = 1.0e-6;

	//smoothness indicators
	double s1 = 13.0 / 12.0 * pow(v1 - 2.0 * v2 + v3, 2) + 0.25 * pow(v1 - 4.0 * v2 + 3.0 * v3, 2);
	double s2 = 13.0 / 12.0 * pow(v2 - 2.0 * v3 + v4, 2) + 0.25 * pow(v2 - v4, 2);
	double s3 = 13.0 / 12.0 * pow(v3 - 2.0 * v4 + v5, 2) + 0.25 * pow(3.0 * v3 - 4.0 * v4 + v5, 2);

	//computing nonlinear weights w1, w2, w3
	double c1 = 0.3 / (pow(eps + s1, 2));
	double c2 = 0.6 / (pow(eps + s2, 2));
	double c3 = 0.1 / (pow(eps + s3, 2));

	double w1 = c1 / (c1 + c2 + c3);
	double w2 = c2 / (c1 + c2 + c3);
	double w3 = c3 / (c1 + c2 + c3);

	//candiate stencils
	double q1 = -v1 / 6.0 + 5.0 / 6.0 * v2 + 1.0 / 3.0 * v3;
	double q2 = v2 / 3.0 + 5.0 / 6.0 * v3 - v4 / 6.0;
	double q3 = 11.0 / 6.0 * v3 - 7.0 / 6.0 * v4 + v5 / 3.0;

	//reconstructed value at interface
	double f = (w1 * q1 + w2 * q2 + w3 * q3);

	return f;
}
