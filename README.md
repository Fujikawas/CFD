# Start here
It's a project mainly with C++ to test typical numerical schemes.

The main function stays in *main.cpp*. The different schemes are written as functions and the whole solving process is in the function completed i.e. nothing in commen is considered and in the main function implemented.

# For solving 1D heat transfer equation $\frac{\partial u}{\partial t}=\alpha \frac{\partial^2 u}{\partial x ^2}$
With the *1D_FTCS.cpp, 1D_BTCS.cpp* and *1D_CrackNicolson.cpp* are the forward Euler, backward Euler and trapezoid discretization 
scheme for time tested.Their declarations are in *1D_BTCS.h*. The central FD scheme is choosed for spatial discretization.

The algebraicEqs.cpp and ~.h is written to solve the algebraic 
equation system with Thomas's algorithm, which is simplified Gauss's Elimination method for EQS with tridiagonal coefficient matrix. They 
are required in the implicit and semi-implicit schemes.

The visulization of the solution and discretization error comparing to the analytical solution is implemented with Python in *readu.py*. 
