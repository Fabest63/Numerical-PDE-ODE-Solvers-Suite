function advection_diffusion_mol_theta
    % Method of Lines (MOL) & Theta-Scheme Solver for 1D Advection-Diffusion PDE
    % PDE: u_t = u_xx - 0.5*u_x + g(x,t) on [0,1]x[0,1]
    % Boundary Conditions: Dirichlet u(0,t) = exp(-t), Neumann u_x(1,t) = 0
    m = input('Enter number of spatial subdivisions m: ');
    n = input('Enter number of temporal subdivisions n: ');
    fprintf('Select time-stepping scheme:\n');
    fprintf('  0 - Explicit Forward Euler (theta = 0.0)\n');
    fprintf('  1 - Crank-Nicolson         (theta = 0.5)\n');
    fprintf('  2 - Implicit Backward Euler(theta = 1.0)\n');
    method_idx = input('Option: ');
    thetas = [0, 0.5, 1];
    theta  = thetas(method_idx + 1);

    h = 1 / m;
    x = linspace(0, 1, m + 1)';
    x_interior = x(2:end);

    tau = 1 / n;
    t = linspace(0, 1, n + 1);

    % Initial condition
    W = initial_condition(x_interior);

    % Assemble sparse finite-difference operators incorporating Neumann boundary at x = 1
    D  = spdiags([ones(m,1), -2*ones(m,1), ones(m,1)], -1:1, m, m);
    A  = spdiags([-ones(m,1), zeros(m,1), ones(m,1)], -1:1, m, m);
    T1 = sparse(m, m); T1(m, m-1) = 1;
    T2 = sparse(m, m); T2(m, m-1) = 1;

    Dh = (1 / h^2) * (D + T1);
    Ah = -(1 / (4 * h)) * (A + T2);
    Ph = sparse(Dh + Ah);

    % Source term and Dirichlet boundary contribution at x = 0
    Bh = zeros(m, n + 1);
    for j = 1:n+1
        Bh(1, j) = ((4 + h) / (4 * h^2)) * exp(-t(j)) + source_term(x_interior(1), t(j));
        for i = 2:m
            Bh(i, j) = source_term(x_interior(i), t(j));
        end
    end

    % Exact solution grid
    exact_grid = zeros(m, n + 1);
    for i = 1:m
        for j = 1:n+1
            exact_grid(i, j) = exact_solution(x_interior(i), t(j));
        end
    end

    % Time-stepping loop via Sparse LU Factorization
    Solu = zeros(m, n + 1);
    Solu(:, 1) = W;
    I = speye(size(Ph));

    system_matrix = I - tau * theta * Ph;
    [L, U] = lu(system_matrix);

    for j = 1:n
        rhs = (I + tau * (1 - theta) * Ph) * W + tau * (1 - theta) * Bh(:, j) + tau * theta * Bh(:, j + 1);
        W = U \ (L \ rhs);
        Solu(:, j + 1) = W;
    end

    % Compute maximum error
    error_abs = max(max(abs(Solu - exact_grid)));
    fprintf('\nMaximum absolute error (m = %d, n = %d, theta = %.1f): %.4e\n', m, n, theta, error_abs);

    % Plot 2D profile at final time t = 1 and export table
    plot_final_profile(x_interior, Solu, exact_grid, m, n, theta);
    export_mol_table(x_interior, m, n, theta, exact_grid, Solu, error_abs);
end

function u = initial_condition(x)
    u = -3 * x.^2 + 6 * x + 1;
end

function val = source_term(x, t)
    val = exp(-t) .* (3 * x.^2 - 9 * x + 8);
end

function u = exact_solution(x, t)
    u = exp(-t) .* (-3 * x.^2 + 6 * x + 1);
end

function plot_final_profile(x_interior, Solu, exact_grid, m, n, theta)
    figure; hold on;
    plot(x_interior, Solu(:, n + 1), 'r.-', 'DisplayName', 'Numerical Approximation');
    plot(x_interior, exact_grid(:, n + 1), 'k--', 'DisplayName', 'Exact Analytical Solution');
    title(sprintf('Advection-Diffusion Solution (m=%d, n=%d, \\theta=%.1f, t=1)', m, n, theta));
    legend('Location', 'best');
    xlabel('Spatial coordinate x');
    ylabel('u(x, 1)');
    grid on;
end

function export_mol_table(x_interior, m, n, theta, exact_grid, Solu, error_abs)
    u = length(x_interior);
    RowNames = {'Exact', 'Approx'};

    fid = fopen('advection_diffusion_results.txt', 'a');
    fprintf(fid, '\n');
    fprintf(fid, 'Spatial subdivisions m:  %d\n', m);
    fprintf(fid, 'Temporal subdivisions n: %d\n', n);
    fprintf(fid, 'Theta parameter:         %.1f\n', theta);

    fprintf(fid, '+------------');
    for k = 1:u; fprintf(fid, '+------------'); end
    fprintf(fid, '+\n');

    fprintf(fid, '| %-10s ', 'x');
    for k = 1:u; fprintf(fid, '| %-10.5f ', x_interior(k)); end
    fprintf(fid, '|\n');

    fprintf(fid, '+------------');
    for k = 1:u; fprintf(fid, '+------------'); end
    fprintf(fid, '+\n');

    Data = [exact_grid(:, n + 1)'; Solu(:, n + 1)'];
    for k = 1:2
        fprintf(fid, '| %-10s ', RowNames{k});
        for p = 1:u; fprintf(fid, '| %-10.3e ', Data(k, p)); end
        fprintf(fid, '|\n');

        fprintf(fid, '+------------');
        for p = 1:u; fprintf(fid, '+------------'); end
        fprintf(fid, '+\n');
    end

    fprintf(fid, 'Maximum Absolute Error: %.4e\n', error_abs);
    fclose(fid);
end
