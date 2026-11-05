# Numerical-PDE-ODE-Solvers-Suite

Numerical analysis implementations developed for the **Numerical Methods in Partial Differential Equations** course (Grade: **9.00/10 — Distinction**) within the BSc in Mathematics at the University of La Laguna.

*(Note: For the nonlinear Fisher-KPP reaction-diffusion equation solved via L-stable Rosenbrock schemes, see my dedicated repository: [`Fisher-KPP-Rosenbrock-Solver`](https://github.com/Fabest63/Fisher-KPP-Rosenbrock-Solver)).*

## Implemented Solvers
- **1D Advection-Diffusion PDE via Method of Lines & $\theta$-Schemes (`practica3_adveccion_difusion.m`, `practica3_pdepe.m`)**:
  Solves $u_t = u_{xx} - \frac{1}{2}u_x + g(x,t)$ under mixed Dirichlet-Neumann boundary conditions. Constructs sparse central-difference operators (`spdiags`) and integrates in time via Explicit Euler ($\theta=0$), **Crank-Nicolson ($\theta=0.5$)**, and Implicit Euler ($\theta=1$) using sparse $LU$ factorization (`[L,U] = lu(A)`), benchmarked against `pdepe` and exact analytical solutions.
- **Finite Element Method (FEM) with Linear & Cubic B-Splines (`pract4.m`, `pract4_B1.m`, `pract4_B3.m`)**:
  Solves the variable-coefficient Sturm-Liouville boundary value problem $-(x^2 u')' + 6u = 6x^3$ on $[1,2]$. Implements $C^2$ cubic B-spline basis functions from scratch with modified boundary splines, assembling a 7-band sparse mass/stiffness matrix via adaptive quadrature (`integral`) across refined meshes ($n \in \{5, 10, 20, 40\}$).
- **Nonlinear Boundary Value Problem via Newton-Raphson Shooting Method (`pract2.m`)**:
  Solves $y'' = -2yy'/t$ by coupling a 2nd-order explicit Runge-Kutta (Heun) integrator with the variational sensitivity ODE to update shooting slopes quadratically via Newton-Raphson.
- **Stiff Dynamical Systems (`pract1.m`)**:
  Compares explicit and implicit time-stepping schemes on stiff coupled RLC differential equations (stiffness ratio $> 1.6 \times 10^4$).
