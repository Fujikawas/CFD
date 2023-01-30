#include "1D_BTCS.h"

using namespace std;

void CRWENO5_Dirichlet()
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

		rhs_CR(nx, dx, u[j - 1], r);
		for (int i = 1; i < nx; i++)
		{
			ut[i] = u[j - 1][i] + dt * r[i];
		}
		rhs_CR(nx, dx, ut, r);
		for (int i = 1; i < nx; i++)
		{
			ut[i] = 0.75 * u[j - 1][i] + 0.25 * ut[i] + dt / 4.0 * r[i];
		}
		rhs_CR(nx, dx, ut, r);
		for (int i = 1; i < nx; i++)
		{
			u[j][i] = u[j - 1][i] / 3.0 + 2.0 / 3.0 * ut[i] + dt * 2.0 / 3.0 * r[i];
		}

	}

	ofstream outfile("CRWENO5.dat");
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

void rhs_CR(int nx, double dx, vector<double> u, vector<double>& r)
{
	vector<double> uL(nx, 0);
	vector<double> uR(nx + 1, 0);

	wenoL_CR(nx, u, uL);
	wenoR_CR(nx, u, uR);

	for (int i = 1; i < nx; i++)
	{
		if (u[i] >= 0.0)
			r[i] = -u[i] * (uL[i] - uL[i - 1]) / dx;
		else
			r[i] = -u[i] * (uR[i + 1] - uR[i]) / dx;
	}
	return;
}

void wenoL_CR(int nx, vector<double> u, vector<double>& uL)
{
	vector<double> a(nx, 0);
	vector<double> b(nx, 0);
	vector<double> c(nx, 0);
	vector<double> r(nx, 0);

	int i = 0;
	b[i] = 2.0 / 3.0;
	c[i] = 1.0 / 3.0;
	r[i] = (u[i] + 5.0 * u[i + 1]) / 6.0;
	
	i = 1;
	double v1 = 2.0 * u[i - 1] - u[i];
	double v2 = u[i - 1];
	double v3 = u[i];
	double v4 = u[i + 1];
	double v5 = u[i + 2];

	vector<double> returnValue_for_a(6,0);
	returnValue_for_a = wcL_CR(v1, v2, v3, v4, v5);
	a[i] = returnValue_for_a[0];
	b[i] = returnValue_for_a[1];
	c[i] = returnValue_for_a[2];
	r[i] = returnValue_for_a[3] * u[i - 1] + returnValue_for_a[4] * u[i] + returnValue_for_a[5] * u[i + 1];

	for (i = 2; i < nx - 1; i++)
	{
		v1 = u[i - 2];
		v2 = u[i - 1];
		v3 = u[i];
		v4 = u[i + 1];
		v5 = u[i + 2];
		returnValue_for_a = wcL_CR(v1, v2, v3, v4, v5);
		a[i] = returnValue_for_a[0];
		b[i] = returnValue_for_a[1];
		c[i] = returnValue_for_a[2];
		r[i] = returnValue_for_a[3] * u[i - 1] + returnValue_for_a[4] * u[i] + returnValue_for_a[5] * u[i + 1];
	}
	i = nx - 1;
	a[i] = 1.0 / 3.0;
	b[i] = 2.0 / 3.0;
	r[i] = (5.0 * u[i] + u[i + 1]) / 6.0;

	thomasTridiagonal_CRWENO(a, b, c, uL, r, 1, nx);
	return;
}
void wenoR_CR(int nx, vector<double> u, vector<double>& uR)
{
	vector<double> a(nx+1, 0);
	vector<double> b(nx+1, 0);
	vector<double> c(nx+1, 0);
	vector<double> r(nx+1, 0);
	vector<double> returnValue_for_a(6, 0);

	int i = 1;
	b[i] = 2.0 / 3.0;
	c[i] = 1.0 / 3.0;
	r[i] = (u[i-1] + 5.0 * u[i]) / 6.0;

	double v1,v2,v3,v4,v5 = 0;

	for (i = 2; i < nx - 1; i++)
	{
		v1 = u[i - 2];
		v2 = u[i - 1];
		v3 = u[i];
		v4 = u[i + 1];
		v5 = u[i + 2];
		returnValue_for_a = wcR_CR(v1, v2, v3, v4, v5);
		a[i] = returnValue_for_a[0];
		b[i] = returnValue_for_a[1];
		c[i] = returnValue_for_a[2];
		r[i] = returnValue_for_a[3] * u[i - 1] + returnValue_for_a[4] * u[i] + returnValue_for_a[5] * u[i + 1];
	}
	i = nx - 1;
	v1 = u[i - 2];
	v2 = u[i - 1];
	v3 = u[i];
	v4 = u[i + 1];
	v5 = 2.0 * u[i + 1] - u[i];

	returnValue_for_a = wcR_CR(v1, v2, v3, v4, v5);
	a[i] = returnValue_for_a[0];
	b[i] = returnValue_for_a[1];
	c[i] = returnValue_for_a[2];
	r[i] = returnValue_for_a[3] * u[i - 1] + returnValue_for_a[4] * u[i] + returnValue_for_a[5] * u[i + 1];

	i = nx;
	a[i] = 1.0 / 3.0;
	b[i] = 2.0 / 3.0;
	r[i] = (5.0 * u[i - 1] + u[i]) / 6.0;

	thomasTridiagonal_CRWENO(a, b, c, uR, r, 2, nx+1);

	return;
}

vector<double> wcL_CR(double v1, double v2, double v3, double v4, double v5)
{
	double eps = 1.0e-6;
	vector<double> returnValue(6, 0);
	//smoothness indicators
	double s1 = 13.0 / 12.0 * pow(v1 - 2.0 * v2 + v3, 2) + 0.25 * pow(v1 - 4.0 * v2 + 3.0 * v3, 2);
	double s2 = 13.0 / 12.0 * pow(v2 - 2.0 * v3 + v4, 2) + 0.25 * pow(v2 - v4, 2);
	double s3 = 13.0 / 12.0 * pow(v3 - 2.0 * v4 + v5, 2) + 0.25 * pow(3.0 * v3 - 4.0 * v4 + v5, 2);

	//computing nonlinear weights w1, w2, w3
	double c1 = 2.0e-1 / (pow(eps + s1, 2));
	double c2 = 5.0e-1 / (pow(eps + s2, 2));
	double c3 = 3.0e-1 / (pow(eps + s3, 2));

	double w1 = c1 / (c1 + c2 + c3);
	double w2 = c2 / (c1 + c2 + c3);
	double w3 = c3 / (c1 + c2 + c3);

	//candiate stencils
	returnValue[0] = (2.0 * w1 + w2) / 3.0;
	returnValue[1] = (w1 + 2.0 * w2 + 2.0 * w3) / 3.0;
	returnValue[2] = w3 / 3.0;

	returnValue[3] = w1 / 6.0;
	returnValue[4] = (5.0 * w1 + 5.0 * w2 + w3) / 6.0;
	returnValue[5] = (w2 + 5.0 * w3) / 6.0;

	return returnValue;
}
vector<double> wcR_CR(double v1, double v2, double v3, double v4, double v5)
{
	double eps = 1.0e-6;
	vector<double> returnValue(6, 0);

	//smoothness indicators
	double s1 = 13.0 / 12.0 * pow(v1 - 2.0 * v2 + v3, 2) + 0.25 * pow(v1 - 4.0 * v2 + 3.0 * v3, 2);
	double s2 = 13.0 / 12.0 * pow(v2 - 2.0 * v3 + v4, 2) + 0.25 * pow(v2 - v4, 2);
	double s3 = 13.0 / 12.0 * pow(v3 - 2.0 * v4 + v5, 2) + 0.25 * pow(3.0 * v3 - 4.0 * v4 + v5, 2);

	//computing nonlinear weights w1, w2, w3
	double c1 = 0.3 / (pow(eps + s1, 2));
	double c2 = 0.5 / (pow(eps + s2, 2));
	double c3 = 0.2 / (pow(eps + s3, 2));

	double w1 = c1 / (c1 + c2 + c3);
	double w2 = c2 / (c1 + c2 + c3);
	double w3 = c3 / (c1 + c2 + c3);

	returnValue[0] = w1 / 3.0;
	returnValue[1] = (w3 + 2.0 * w2 + 2.0 * w1) / 3.0;
	returnValue[2] = (2.0 * w3 + w2) / 3.0;

	returnValue[3] = (w2 + 5.0 * w1) / 6.0;
	returnValue[4] = (5.0 * w3 + 5.0 * w2 + w1) / 6.0;
	returnValue[5] = w3 / 6.0;

	return returnValue;
}
void thomasTridiagonal_CRWENO(vector<double> a, vector<double> b, vector<double> c, vector<double>& u, vector<double> q, int s, int e)
{
	vector<double> m(e, 0);
	u[s-1] = q[s-1] / b[s-1];

	for (int i = s; i < e; i++)
	{
		m[i] = c[i - 1] / b[i - 1];
		b[i] = b[i] - m[i] * a[i];
		u[i] = (q[i] - u[i - 1] * a[i]) / b[i];
	}
	u[e - 1] = q[e - 1];
	for (int i = e - 2; i >= 0; i--)
	{
		u[i] = u[i] - u[i + 1] * m[i + 1];
	}
	return;
}