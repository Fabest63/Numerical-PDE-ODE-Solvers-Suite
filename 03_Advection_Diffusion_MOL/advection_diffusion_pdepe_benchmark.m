function advection_diffusion_pdepe_benchmark
    % 3D Benchmark of the 1D Advection-Diffusion Problem using MATLAB's pdepe solver
    m = input('Enter number of spatial subdivisions m: ');
    n = input('Enter number of temporal subdivisions n: ');

    x = linspace(0, 1, m + 1);
    t = linspace(0, 1, n + 1);

    sol = pdepe(0, @pde_equation, @initial_condition, @boundary_conditions, x, t);

    % Exact analytical solution over space-time mesh
    [X, T] = meshgrid(x, t);
    U_exact = exp(-T) .* (-3 * X.^2 + 6 * X + 1);

    % Maximum error across the space-time grid
    max_error = max(max(abs(sol - U_exact)));
    fprintf('Maximum space-time error using pdepe: %.4e\n', max_error);

    % 3D Surface Plots
    figure;
    surf(x, t, sol);
    title('Numerical Solution Surface u(x,t) via pdepe');
    xlabel('Space x'); ylabel('Time t'); zlabel('u(x,t)');
    shading interp; colorbar;

    figure;
    surf(x, t, abs(sol - U_exact));
    title('Absolute Error Surface |u_{num}(x,t) - u_{exact}(x,t)|');
    xlabel('Space x'); ylabel('Time t'); zlabel('Absolute Error');
    shading interp; colorbar;
end

function [c, f, s] = pde_equation(x, t, ~, dudx)
    c = 1;
    f = dudx;
    s = -0.5 * dudx + source_term(x, t);
end

function u = initial_condition(x)
    u = -3 * x.^2 + 6 * x + 1;
end

function [pl, ql, pr, qr] = boundary_conditions(~, ul, ~, ~, t)
    % Left boundary: Dirichlet u(0,t) = exp(-t)
    pl = ul - exp(-t);
    ql = 0;

    % Right boundary: Homogeneous Neumann u_x(1,t) = 0
    pr = 0;
    qr = 1;
end

function val = source_term(x, t)
    val = exp(-t) .* (3 * x.^2 - 9 * x + 8);
end
