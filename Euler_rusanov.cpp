#include "1D_BTCS.h"

using namespace std;

void Euler_rusanov()
{
	int nx = 256;
	double dx = 1.0 / nx;
	double dt = 0.0001;
	double t = 0.2;
	int nt = ceil(t / dt);
	vector<double> x(nx, 0);
	vector<vector<double>> qt(3, vector<double>(nx, 0));//temperory array by RK3 scheme
	vector<vector<vector<double>>> q(nt, vector<vector<double>>(3, vector<double>(nx, 0)));//all timesteps
	vector<vector<double>> r(3, vector<double>(nx, 0));//general spatial FD
	double gamma = 1.4;//specific gas ratio

	//Sod's Riemann problem
	//left side
	double rhoL = 1.0;
	double uL = 0.0;
	double pL = 1.0;
	//right side
	double rhoR = 0.125;
	double uR = 0.0;
	double pR = 0.1;

	double rho, u, p, e = 0;
	//grid
	for (int i = 0; i < nx; i++)
	{
		x[i] = (i + 0.5) * dx;
	}
	double xc = 0.5;
	for (int i = 0; i < nx; i++)
	{
		rho = x[i] > 0.5 ? rhoR : rhoL;
		u = x[i] > 0.5 ? uR : uL;
		p = x[i] > 0.5 ? pR : pL;
		e = p / (rho * (gamma - 1.0)) + 0.5 * u * u;
		q[0][0][i] = rho;
		q[0][1][i] = rho * u;
		q[0][2][i] = rho * e;
	}

	for (int j = 1; j < nt; j++)// one time step
	{
		rhs_rusanov(nx, dx, gamma, q[j - 1], r);
		for (int i = 0; i < nx; i++)
		{
			qt[0][i] = q[j - 1][0][i] + dt * r[0][i];
			qt[1][i] = q[j - 1][1][i] + dt * r[1][i];
			qt[2][i] = q[j - 1][2][i] + dt * r[2][i];
		}
		rhs_rusanov(nx, dx, gamma, qt, r);
		for (int i = 0; i < nx; i++)
		{
			qt[0][i] = 0.75 * q[j - 1][0][i] + 0.25 * qt[0][i] + dt / 4.0 * r[0][i];
			qt[1][i] = 0.75 * q[j - 1][1][i] + 0.25 * qt[1][i] + dt / 4.0 * r[1][i];
			qt[2][i] = 0.75 * q[j - 1][2][i] + 0.25 * qt[2][i] + dt / 4.0 * r[2][i];
		}
		rhs_rusanov(nx, dx, gamma, qt, r);
		for (int i = 0; i < nx; i++)
		{
			q[j][0][i] = q[j - 1][0][i] / 3.0 + 2.0 / 3.0 * qt[0][i] + dt * 2.0 / 3.0 * r[0][i];
			q[j][1][i] = q[j - 1][1][i] / 3.0 + 2.0 / 3.0 * qt[1][i] + dt * 2.0 / 3.0 * r[1][i];
			q[j][2][i] = q[j - 1][2][i] / 3.0 + 2.0 / 3.0 * qt[2][i] + dt * 2.0 / 3.0 * r[2][i];
		}
	}

	ofstream outfile("rusanov.dat");
	if (outfile.is_open())
	{
		for (int m = 0; m < 3; m++)
		{
			for (int i = 0; i < nx; i++)
			{
				outfile << q[nt - 1][m][i] << " ";
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

void rhs_rusanov(int nx, double dx, double gamma, vector<vector<double>> q, vector<vector<double>>& r)
{
	//left and right side fluxes at the interface
	vector<vector<double>> qL(3, vector<double>(nx + 1, 0));
	vector<vector<double>> qR(3, vector<double>(nx + 1, 0));

	vector<vector<double>> fL(3, vector<double>(nx + 1, 0));
	vector<vector<double>> fR(3, vector<double>(nx + 1, 0));
	vector<vector<double>> f(3, vector<double>(nx + 1, 0));

	qL = wenoL_roe(nx, q);
	qR = wenoR_roe(nx, q);

	flux_roe(nx, gamma, qL, fL);
	flux_roe(nx, gamma, qR, fR);

	rusanov(nx, gamma, qL, qR, f, fL, fR);
	for (int i = 0; i < nx; i++)
	{
		for (int m = 0; m < 3; m++)
		{
			r[m][i] = -(f[m][i + 1] - f[m][i]) / dx;
		}
	}
	return;
}
void rusanov(int nx, double gamma, vector<vector<double>> qL, vector<vector<double>> qR, vector<vector<double>>& f, vector<vector<double>> fL, vector<vector<double>> fR)
{

	vector<double> ps(nx + 1, 0);
	vector<double> rad(nx, 0);
	//spectral radius of Jacobian
	double gm = gamma - 1.0;
	for (int i = 0; i < nx + 1; i++)
	{
		//Leftand right states :
		double rhLL = qL[0][i];
		double uuLL = qL[1][i] / rhLL;
		double eeLL = qL[2][i] / rhLL;
		double ppLL = gm * (eeLL * rhLL - 0.5 * rhLL * (uuLL * uuLL));
		double hhLL = eeLL + ppLL / rhLL;


		double rhRR = qR[0][i];
		double uuRR = qR[1][i] / rhRR;
		double eeRR = qR[2][i] / rhRR;
		double ppRR = gm * (eeRR * rhRR - 0.5 * rhRR * (uuRR * uuRR));
		double hhRR = eeRR + ppRR / rhRR;

		double alpha = 1.0 / (sqrt(abs(rhLL)) + sqrt(abs(rhRR)));

		double uu = (sqrt(abs(rhLL)) * uuLL + sqrt(abs(rhRR)) * uuRR) * alpha;
		double hh = (sqrt(abs(rhLL)) * hhLL + sqrt(abs(rhRR)) * hhRR) * alpha;
		double aa = sqrt(abs(gm * (hh - 0.5 * uu * uu)));

		ps[i] = abs(aa + uu);
	}			

	for (int i = 0; i < nx + 1; i++)
	{
			for (int m = 0; m < 3; m++)
			{
				f[m][i] = 0.5 * (fR[m][i] + fL[m][i]) - 0.5 * ps[i] * (qR[m][i] - qL[m][i]); 
			}
		
	}
}

