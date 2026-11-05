function fem_sturm_liouville_bspline_main
    % Main Launcher: Finite Element Method (Rayleigh-Ritz) with B-Spline Basis Functions
    % Solves -(x^2 * u')' + 6u = 6x^3 on [1, 2] with homogeneous Dirichlet BCs u(1) = u(2) = 0
    order = input('Select B-Spline basis order (1 for Linear B1, 3 for Cubic B3): ');
    if order == 1
        fem_bspline_linear;
    elseif order == 3
        fem_bspline_cubic;
    else
        fprintf('Invalid selection. Please enter 1 (Linear) or 3 (Cubic).\n');
    end
end

