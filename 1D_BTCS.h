#pragma once
constexpr auto Pi = 3.1415926;

#include <iostream>
#include <fstream>
#include <ostream>
#include <iterator>
#include <vector>

using namespace std;

//for heat transfer
void BTCS_1D();
void FTCS_1D();
void CNTCS_1D();
void RK2_1D();
void RK3_1D();
void RK4_1D();

//for solution of algebraic eqs.
void thomasTridiagonal(vector<double>, vector<double>, vector<double>, vector<double>&, vector<double>);

//Burger's equation
// for WENO5
void WENO5_Dirichlet();
void rhs(int, double, vector<double>, vector<double>&);
void wenoL(int, vector<double>, vector<double>&);
void wenoR(int, vector<double>, vector<double>&);
double wcL(double, double, double, double, double);
double wcR(double, double, double, double, double);
//for CRWENO5
void CRWENO5_Dirichlet();
void rhs_CR(int, double, vector<double>, vector<double>&);// Calculate right hand term of the inviscid Burgers equation
void wenoL_CR(int, vector<double>, vector<double>&);// CRWENO reconstruction ofr upwind direction (positive and left to right)
void wenoR_CR(int, vector<double>, vector<double>&);
vector<double> wcL_CR(double, double, double, double, double);//nonlinear weights for upwind direction
vector<double> wcR_CR(double, double, double, double, double);
void thomasTridiagonal_CRWENO(vector<double>, vector<double>, vector<double>, vector<double>&, vector<double>,int, int);
// for FVM-Flux splitting
void FluxSplitting_Burgers();
void rhs_FS(int, double, vector<double>, vector<double>&);
void waveSpeed(int, vector<double>, vector<double>&);
vector<double> wenoL_FS(int, vector<double>);
vector<double> wenoR_FS(int, vector<double>);
//for FVM-Riemann solver
void Riemann_period();
void rhs_RM(int, double, vector<double>, vector<double>&);
void Riemann(int, vector<double>, vector<double>, vector<double>, vector<double>&, vector<double>, vector<double>);

//1D Euler Solver
void Euler_roe();
void rhs_roe(int, double, double, vector<vector<double>>, vector<vector<double>>&);
void flux_roe(int, double, vector<vector<double>>, vector<vector<double>>&);
void roe(int, double, vector<vector<double>>, vector<vector<double>>, vector<vector<double>>&,vector<vector<double>>, vector<vector<double>>);
vector<vector<double>> wenoL_roe(int, vector<vector<double>>);
vector<vector<double>> wenoR_roe(int, vector<vector<double>>);
//HLLC
void Euler_hllc();
void rhs_hllc(int, double, double, vector<vector<double>>, vector<vector<double>>&);
void hllc(int, double, vector<vector<double>>, vector<vector<double>>, vector<vector<double>>&, vector<vector<double>>, vector<vector<double>>);
//Rusanov
//HLLC
void Euler_rusanov();
void rhs_rusanov(int, double, double, vector<vector<double>>, vector<vector<double>>&);
void rusanov(int, double, vector<vector<double>>, vector<vector<double>>, vector<vector<double>>&, vector<vector<double>>, vector<vector<double>>);