# Start here
It's a project mainly with C++ to test typical numerical schemes. The most interpretations to the algorithms in this project could be found in: Pawar, S.; San, O. CFD Julia: A Learning Module Structuring an Introductory Course on Computational Fluid Dynamics. Fluids 2019, 4, 159. https://doi.org/10.3390/fluids4030159

The main function stays in *main.cpp*. The different schemes are written as functions and the whole solving process is in the function completed i.e. nothing in commen is considered and in the main function implemented.

# For solving 1D heat transfer equation $\frac{\partial u}{\partial t}=\alpha \frac{\partial^2 u}{\partial x ^2}$
With the *1D_FTCS.cpp, 1D_BTCS.cpp* and *1D_CrackNicolson.cpp* are the forward Euler, backward Euler and trapezoid discretization 
scheme for time tested. The three typical Runge-Kutta schemes(in three orders) for time dimension are also tested. Their declarations are in *1D_BTCS.h*. The central FD scheme is choosed for spatial discretization.

The algebraicEqs.cpp and ~.h is written to solve the algebraic 
equation system with Thomas's algorithm, which is simplified Gauss's Elimination method for EQS with tridiagonal coefficient matrix. They 
are required in the implicit and semi-implicit schemes.

The visualization of the solution and discretization error comparing to the analytical solution is implemented with Python in *readu.py*. 
# For solving 1D inviscid Burgers equation $\frac{\partial u}{\partial t}+u\frac{\partial u}{\partial x}=0$
The source codes with input, algorithm and output are in Burgers_*.cpp. The WENO5 and CRWENO5 are FDM with WENO and CRWENO scheme in 5th order and the boundary condition are Dirichlet type. 

The FluxSplitting and Riemeann are FVM with periodic condition at the boundary. And the first one use Lax-Friedrichs Flux spltting method to calculate the positive and negative flux components at a cell center, while the second one use Rusanov scheme as the Riemann solver to compute the flux on the cell surfaces.
# For solving 1D Euler equation $\frac{\partial q}{\partial t}+\frac{\partial F}{\partial x}=0$ with $q=(\rho\quad \rho u\quad \rho e)^T$ and $F=(\rho u\quad \rho u^2 +p\quad \rho uh)^T$ in which $h=e+\frac{p}{\rho}$, $p=\rho(\gamma-1)(e-\frac{1}{2}u^2)$
The source codes with input, algorithm and output are in Euler_*.cpp. 

The third-order Runge-Kutta numerical method is used for time integration and the WENO-5 reconstruction is used to compute the left and sight states of the fluxes at the interface.

Three different Riemann solvers to compute the flux on the cell faces.


# For solving 2D poisson equation $\frac{\partial^2 u}{\partial x ^2}+\frac{\partial^2 u}{\partial y ^2}=f$
We perform the assessment of direct solver using the method of manufactured solution. We assume certain field u and compute the source term f at each grid location. We then solve the Poisson equation for this source term f and compare the numerically calculated field u with the exact solution field u.

A direct Poisson solver is implemented with FFT with FFTW module(Poisson_FFT.cpp). Note that you need to install the FFTW module at first to run the code.

The iterative solver are with stationary method like Gauss-Seidel method and non-stationary method like CG method.

Furthermore, a multigrid solver(two layers) is implemented.
