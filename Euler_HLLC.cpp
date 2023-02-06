#include "1D_BTCS.h"

using namespace std;

void Euler_hllc()
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
		rhs_hllc(nx, dx, gamma, q[j - 1], r);
		for (int i = 0; i < nx; i++)
		{
			qt[0][i] = q[j - 1][0][i] + dt * r[0][i];
			qt[1][i] = q[j - 1][1][i] + dt * r[1][i];
			qt[2][i] = q[j - 1][2][i] + dt * r[2][i];
		}
		rhs_hllc(nx, dx, gamma, qt, r);
		for (int i = 0; i < nx; i++)
		{
			qt[0][i] = 0.75 * q[j - 1][0][i] + 0.25 * qt[0][i] + dt / 4.0 * r[0][i];
			qt[1][i] = 0.75 * q[j - 1][1][i] + 0.25 * qt[1][i] + dt / 4.0 * r[1][i];
			qt[2][i] = 0.75 * q[j - 1][2][i] + 0.25 * qt[2][i] + dt / 4.0 * r[2][i];
		}
		rhs_hllc(nx, dx, gamma, qt, r);
		for (int i = 0; i < nx; i++)
		{
			q[j][0][i] = q[j - 1][0][i] / 3.0 + 2.0 / 3.0 * qt[0][i] + dt * 2.0 / 3.0 * r[0][i];
			q[j][1][i] = q[j - 1][1][i] / 3.0 + 2.0 / 3.0 * qt[1][i] + dt * 2.0 / 3.0 * r[1][i];
			q[j][2][i] = q[j - 1][2][i] / 3.0 + 2.0 / 3.0 * qt[2][i] + dt * 2.0 / 3.0 * r[2][i];
		}
	}

	ofstream outfile("hllc.dat");
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

void rhs_hllc(int nx, double dx, double gamma, vector<vector<double>> q, vector<vector<double>>& r)
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

	hllc(nx, gamma, qL, qR, f, fL, fR);
	for (int i = 0; i < nx; i++)
	{
		for (int m = 0; m < 3; m++)
		{
			r[m][i] = -(f[m][i + 1] - f[m][i]) / dx;
		}
	}
	return;
}
void hllc(int nx, double gamma, vector<vector<double>> uL, vector<vector<double>> uR, vector<vector<double>>& f, vector<vector<double>> fL, vector<vector<double>> fR)
{
	vector<double> Ds(3, 0);
	Ds[1] = 1.0;
	double gm = gamma - 1.0;

	for (int i = 0; i < nx + 1; i++)
	{
		//Leftand right states :
		double rhLL = uL[0][i];
		double uuLL = uL[1][i] / rhLL;
		double eeLL = uL[2][i] / rhLL;
		double ppLL = gm * (eeLL * rhLL - 0.5 * rhLL * (uuLL * uuLL));
		double aaLL = sqrt(abs(gamma * ppLL / rhLL));

		double rhRR = uR[0][i];
		double uuRR = uR[1][i] / rhRR;
		double eeRR = uR[2][i] / rhRR;
		double ppRR = gm * (eeRR * rhRR - 0.5 * rhRR * (uuRR * uuRR));
		double aaRR = sqrt(abs(gamma * ppRR / rhRR));

		// compute SLand Sr
		double SL = min(uuLL, uuRR) - max(aaLL, aaRR);
		double SR = max(uuLL, uuRR) + max(aaLL, aaRR);

		// compute compound speed
		double	SP = (ppRR - ppLL + rhLL * uuLL * (SL - uuLL) - rhRR * uuRR * (SR - uuRR)) / (rhLL * (SL - uuLL) - rhRR * (SR - uuRR)); //never get zero

		// compute compound pressure
		double PLR = 0.5 * (ppLL + ppRR + rhLL * (SL - uuLL) * (SP - uuLL) + rhRR * (SR - uuRR) * (SP - uuRR));

		// compute D
		Ds[2] = SP;

		if (SL >= 0.0)
		{
			for (int m = 0; m < 3; m++)
			{
				f[m][i] = fL[m][i];
			}
		}
		else if (SR <= 0.0)
		{
			for (int m = 0; m < 3; m++)
			{
				f[m][i] = fR[m][i];
			}
		}
		else if ((SP >= 0.0) & (SL <= 0.0))
		{
			for (int m = 0; m < 3; m++)
			{
				f[m][i] = (SP * (SL * uL[m][i] - fL[m][i]) + SL * PLR * Ds[m]) / (SL - SP);
			}
		}
		else if ((SP <= 0.0) & (SR >= 0.0))
		{
			for (int m = 0; m < 3; m++)
			{
				f[m][i] = (SP * (SR * uR[m][i] - fR[m][i]) + SR * PLR * Ds[m]) / (SR - SP);
			}
		}
	}
}

