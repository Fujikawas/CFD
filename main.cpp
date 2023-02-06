#include "1D_BTCS.h"

using namespace std;

int main()
{
	//FTCS_1D();
	//BTCS_1D();
	//CNTCS_1D();
	//RK2_1D();
	//RK3_1D();
	//RK4_1D();
	//WENO5_Dirichlet();
	//CRWENO5_Dirichlet();
	//FluxSplitting_Burgers();
	//Riemann_period();
	//Euler_roe();
	//Euler_hllc();
	Euler_rusanov();
	return 0;
}