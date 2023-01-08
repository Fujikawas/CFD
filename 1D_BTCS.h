#pragma once
constexpr auto Pi = 3.1415926;

#include <iostream>
#include <fstream>
#include <ostream>
#include <iterator>
#include <vector>
#include "algebraicEqs.h";

using namespace std;

void BTCS_1D();
void FTCS_1D();
void CNTCS_1D();
void RK2_1D();
void RK3_1D();
void RK4_1D();