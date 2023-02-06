#include "1D_BTCS.h"

using namespace std;

void Euler_roe()
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
		rho = x[i] > 0.5? rhoR:rhoL;
		u = x[i] > 0.5? uR : uL; 
		p = x[i] > 0.5? pR : pL;
		e = p / (rho * (gamma - 1.0)) + 0.5 * u * u;
		q[0][0][i] = rho;
		q[0][1][i] = rho * u;
		q[0][2][i] = rho * e;
	}

	for (int j = 1; j < nt; j++)// one time step
	{
		rhs_roe(nx, dx, gamma, q[j - 1], r);
		for (int i = 0; i < nx; i++)
		{
			qt[0][i] =q[j - 1][0][i] + dt * r[0][i];
			qt[1][i] =q[j - 1][1][i] + dt * r[1][i];
			qt[2][i] =q[j - 1][2][i] + dt * r[2][i];
		}
		rhs_roe(nx, dx, gamma, qt, r);
		for (int i = 0; i < nx; i++)
		{
			qt[0][i] = 0.75 * q[j - 1][0][i] + 0.25 * qt[0][i] + dt / 4.0 * r[0][i];
			qt[1][i] = 0.75 * q[j - 1][1][i] + 0.25 * qt[1][i] + dt / 4.0 * r[1][i];
			qt[2][i] = 0.75 * q[j - 1][2][i] + 0.25 * qt[2][i] + dt / 4.0 * r[2][i];
		}
		rhs_roe(nx, dx, gamma, qt, r);
		for (int i = 0; i < nx; i++)
		{
			q[j][0][i] = q[j - 1][0][i] / 3.0 + 2.0 / 3.0 * qt[0][i] + dt * 2.0 / 3.0 * r[0][i];
			q[j][1][i] = q[j - 1][1][i] / 3.0 + 2.0 / 3.0 * qt[1][i] + dt * 2.0 / 3.0 * r[1][i];
			q[j][2][i] = q[j - 1][2][i] / 3.0 + 2.0 / 3.0 * qt[2][i] + dt * 2.0 / 3.0 * r[2][i];
		}
	}

	ofstream outfile("roe.dat");
	if (outfile.is_open())
	{ 
		for (int m = 0; m < 3; m++)
		{
			for (int i = 0; i < nx; i++)
			{
				outfile << q[nt-1][m][i] << " ";
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

void rhs_roe(int nx, double dx, double gamma, vector<vector<double>> q, vector<vector<double>>& r)
{
	//left and right side fluxes at the interface
	vector<vector<double>> qL(3, vector<double>(nx + 1, 0));
	vector<vector<double>> qR(3, vector<double>(nx + 1, 0));

	vector<vector<double>> fL(3, vector<double>(nx + 1, 0));
	vector<vector<double>> fR(3, vector<double>(nx + 1, 0));
	vector<vector<double>> f(3, vector<double>(nx + 1, 0));

	/*double a = 0;

	vector<double> alp(nx, 0);
	vector<double> alpha(nx, 0);
	for (int i = 0; i < nx; i++)
	{
		double p = (gamma - 1.0) * (q[2][i] - 0.5 * q[1][i] * q[1][i] / q[0][i]);
		f[0][i] = q[1][i];
		f[1][i] = q[1][i] * q[1][i] / q[0][i] + p;
		f[2][i] = q[1][i] * q[2][i] / q[0][i] + p * q[1][i] / q[0][i];
	}
	for (int i = 0; i < nx; i++)
	{
		a = sqrt(gamma * ((gamma - 1) * (q[2][i] - 0.5 * q[1][i] * q[1][i] / q[0][i])) / q[0][i]);
		alp[i] = abs(q[1][i] / q[0][i]) + a;
	}

	for (int i = 2; i < nx - 3; i++)
	{
		alpha[i] = max(alp[i - 2], max(alp[i - 1], max(alp[i], max(alp[i + 1], alp[i + 2]))));
	}
	alpha[0] = alp[0];
	alpha[1] = alp[1];
	alpha[nx - 2] = alp[nx - 2];
	alpha[nx - 1] = alp[nx - 1];*/

	qL = wenoL_roe(nx, q);
	qR = wenoR_roe(nx, q);

	flux_roe(nx, gamma, qL, fL);
	flux_roe(nx, gamma, qR, fR);

	roe(nx, gamma, qL, qR, f, fL,fR);
	for (int i = 0; i < nx; i++)
	{
		for (int m = 0; m < 3; m++)
		{
			r[m][i] = -(f[m][i + 1] - f[m][i]) / dx;
		}
	}
	return;
}
void flux_roe(int nx, double gamma, vector<vector<double>> q, vector<vector<double>>& f)
{
	for (int i = 0; i < nx+1; i++)
	{
		double p = (gamma - 1.0) * (q[2][i] - 0.5 * q[1][i] * q[1][i] / q[0][i]);
		f[0][i] = q[1][i];
		f[1][i] = q[1][i] * q[1][i] / q[0][i] + p;
		f[2][i] = q[1][i] * q[2][i] / q[0][i] + p* q[1][i] / q[0][i];
	}
	return;
}
void roe(int nx, double gamma, vector<vector<double>> uL, vector<vector<double>> uR, vector<vector<double>>& f, vector<vector<double>> fL, vector<vector<double>> fR)
{
	vector<double> dd(3, 0); 
	vector<double> dF(3, 0); 
	vector<double> v(3, 0); 
	double gm = gamma - 1.0;

	for (int i = 0; i < nx+1; i++)
	{
	//Leftand right states :
		double rhLL = uL[0][i];
		double uuLL = uL[1][i] / rhLL;
		double eeLL = uL[2][i] / rhLL;
		double ppLL = gm * (eeLL * rhLL - 0.5 * rhLL * (uuLL * uuLL));
		double hhLL = eeLL + ppLL / rhLL;

		double rhRR = uR[0][i];
		double uuRR = uR[1][i] / rhRR;
		double eeRR = uR[2][i] / rhRR;
		double ppRR = gm * (eeRR * rhRR - 0.5 * rhRR * (uuRR * uuRR));
		double hhRR = eeRR + ppRR / rhRR;

		double alpha = 1.0 / (sqrt(abs(rhLL)) + sqrt(abs(rhRR)));

		double uu = (sqrt(abs(rhLL)) * uuLL + sqrt(abs(rhRR)) * uuRR) * alpha;
		double hh = (sqrt(abs(rhLL)) * hhLL + sqrt(abs(rhRR)) * hhRR) * alpha;
		double aa = sqrt(abs(gm * (hh - 0.5 * uu * uu)));

		double D11 = abs(uu);
		double D22 = abs(uu + aa);
		double D33 = abs(uu - aa);

		double beta = 0.5 / (aa * aa);
		double phi2 = 0.5 * gm * uu * uu;

		//Right eigenvector matrix
		double R11 = 1.0;
		double R21 = uu;
		double R31 = phi2/gm;
		double R12 = beta;
		double R22 = beta * (uu + aa);
		double R32 = beta * (hh + uu * aa);
		double R13 = beta;
		double R23 = beta * (uu - aa);
		double R33 = beta * (hh - uu * aa);

		//Left eigenvector matrix
		double L11 = 1.0 - phi2 / (aa * aa);
		double L12 = gm * uu / (aa * aa);
		double L13 = -gm / (aa * aa);
		double L21= phi2 - uu * aa;
		double L22 = aa - gm * uu;
		double L23 = gm;
		double L31 = phi2 + uu * aa;
		double L32 = -aa - gm * uu;
		double L33 = gm;

		for (int m = 0; m < 3; m++)
		{
			v[m] = 0.5 * (uR[m][i] - uL[m][i]);
		}
		dd[0] = D11 * (L11 * v[0] + L12 * v[1] + L13 * v[2]);
		dd[1] = D22 * (L21 * v[0] + L22 * v[1] + L23 * v[2]);
		dd[2] = D33 * (L31 * v[0] + L32 * v[1] + L33 * v[2]);

		dF[0] = R11 * dd[0] + R12 * dd[1] + R13 * dd[2];
		dF[1] = R21 * dd[0] + R22 * dd[1] + R23 * dd[2];
		dF[2] = R31 * dd[0] + R32 * dd[1] + R33 * dd[2];
		for (int m = 0; m < 3; m++)
		{
			f[m][i] = 0.5 * (fR[m][i] + fL[m][i]) - dF[m];
		}
	}
}

vector<vector<double>> wenoL_roe(int nx, vector<vector<double>> q)
{
	vector<vector<double>> f(3, vector<double>(nx + 1,0));
	for (int m = 0; m < 3; m++)
	{
		int i = -1;
		double	v1 = q[m][i + 3];
		double 	v2 = q[m][i + 2];
		double	v3 = q[m][i + 1];
		double	v4 = q[m][i + 1];
		double 	v5 = q[m][i + 2];
		f[m][i + 1] = wcL(v1, v2, v3, v4, v5);

		i = 0;
		v1 = q[m][i + 1];
		v2 = q[m][i];
		v3 = q[m][i];
		v4 = q[m][i + 1];
		v5 = q[m][i + 2];
		f[m][i + 1] = wcL(v1, v2, v3, v4, v5);

		i = 1;
		v1 = q[m][i - 1];
		v2 = q[m][i - 1];
		v3 = q[m][i];
		v4 = q[m][i + 1];
		v5 = q[m][i + 2];
		f[m][i + 1] = wcL(v1, v2, v3, v4, v5);

		for (i = 2; i < nx - 2; i++)
		{
			v1 = q[m][i - 2];
			v2 = q[m][i - 1];
			v3 = q[m][i];
			v4 = q[m][i + 1];
			v5 = q[m][i + 2];
			f[m][i + 1] = wcL(v1, v2, v3, v4, v5);
		}
		i = nx - 2;
		v1 = q[m][i - 2];
		v2 = q[m][i - 1];
		v3 = q[m][i];
		v4 = q[m][i + 1];
		v5 = q[m][i + 1];
		f[m][i + 1] = wcL(v1, v2, v3, v4, v5);

		i = nx - 1;
		v1 = q[m][i - 2];
		v2 = q[m][i - 1];
		v3 = q[m][i];
		v4 = q[m][i];
		v5 = q[m][i - 1];
		f[m][i + 1] = wcL(v1, v2, v3, v4, v5);
	}
	return f;
}
vector<vector<double>> wenoR_roe(int nx, vector<vector<double>> q)
{
	vector<vector<double>> f(3, vector<double>(nx + 1, 0));
	for (int m = 0; m < 3; m++)
	{
		int i = 0;
		double v1 = q[m][i + 1];
		double v2 = q[m][i];
		double v3 = q[m][i];
		double v4 = q[m][i + 1];
		double v5 = q[m][i + 2];
		f[m][i] = wcR(v1, v2, v3, v4, v5);

		i = 1;
		v1 = q[m][i - 1];
		v2 = q[m][i - 1];
		v3 = q[m][i];
		v4 = q[m][i + 1];
		v5 = q[m][i + 2];
		f[m][i] = wcR(v1, v2, v3, v4, v5);

		for (i = 2; i < nx - 2; i++)
		{
			v1 = q[m][i - 2];
			v2 = q[m][i - 1];
			v3 = q[m][i];
			v4 = q[m][i + 1];
			v5 = q[m][i + 2];
			f[m][i] = wcR(v1, v2, v3, v4, v5);
		}


		i = nx - 2;
		v1 = q[m][i - 2];
		v2 = q[m][i - 1];
		v3 = q[m][i];
		v4 = q[m][i + 1];
		v5 = q[m][i + 1];
		f[m][i] = wcR(v1, v2, v3, v4, v5);

		i = nx - 1;
		v1 = q[m][i - 2];
		v2 = q[m][i - 1];
		v3 = q[m][i];
		v4 = q[m][i];
		v5 = q[m][i - 1];
		f[m][i] = wcR(v1, v2, v3, v4, v5);

		i = nx;
		v1 = q[m][i - 2];
		v2 = q[m][i - 1];
		v3 = q[m][i - 1];
		v4 = q[m][i - 2];
		v5 = q[m][i - 3];
		f[m][i] = wcR(v1, v2, v3, v4, v5);
	}
	return f;
}