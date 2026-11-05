function nonlinear_bvp_shooting_newton
    % Nonlinear Boundary Value Problem (BVP) Solver via Newton-Raphson Shooting Method
    % BVP: y'' = -2*y*y'/t on [1, 2], with y(1) = 1/2, y(2) = 2/3
    l    = input('Enter grid refinement level k (integer between 2 and 6, where n = 2^k): ');
    s0   = input('Enter initial shooting slope guess s0: ');
    nmax = input('Enter maximum number of Newton iterations: ');
    tol  = input('Enter convergence tolerance (e.g., 1e-8): ');

    valid_k = [2, 3, 4, 5, 6];
    if ~ismember(l, valid_k)
        disp('Invalid value for k. Please choose an integer between 2 and 6.');
        return;
    end

    alpha = 1/2;
    beta  = 2/3;
    a = 1;
    b = 2;
    n = 2^l;
    tau = (b - a) / n;
    t = linspace(a, b, n + 1);
    s = s0;
    i = 0;
    Iter = zeros(nmax, n + 1);

    while i <= nmax - 1
        i = i + 1;

        % 1. Solve main Initial Value Problem IVP(s) via 2nd-order Runge-Kutta (Heun)
        y0 = [alpha; s];
        y  = y0;
        Y  = zeros(2, n + 1);
        Y(:, 1) = y0;
        for k = 0:n-1
            f1   = ode_rhs(t(k + 1), y);
            y1   = y + tau * f1;
            f2   = ode_rhs(t(k + 2), y1);
            ytau = y + tau * (f1 + f2) / 2;
            Y(:, k + 2) = ytau;
            y = ytau;
        end

        % 2. Solve Variational Sensitivity IVP(xi) for Newton-Raphson update
        z0 = [0; 1];
        z  = z0;
        Z  = zeros(2, n + 1);
        Z(:, 1) = z0;
        for j = 0:n-1
            g1   = variational_rhs(t(j + 1), z, Y(:, j + 1));
            z1   = z + tau * g1;
            g2   = variational_rhs(t(j + 2), z1, Y(:, j + 2));
            ztau = z + tau * (g1 + g2) / 2;
            Z(:, j + 2) = ztau;
            z = ztau;
        end

        % Newton-Raphson step
        delta_s = -(Y(1, n + 1) - beta) / z(1);
        Iter(i, :) = Y(1, :);
        Iter(all(abs(Iter) < tol, 2), :) = [];

        % Stopping criterion
        if abs(delta_s) <= tol
            Iter(end, :) = [];
            sol    = Y(1, :);
            exact  = t ./ (1 + t);
            error  = max(abs(sol - exact));

            disp(['Exact solution at grid nodes:        ', num2str(exact)]);
            disp(['Numerical solution at grid nodes:    ', num2str(sol)]);
            disp(['Maximum absolute error:              ', num2str(error)]);
            disp(['Total Newton shooting iterations:    ', num2str(i)]);

            plot_shooting_results(i, t, sol, exact, Iter);
            export_shooting_table(t, sol, exact, error, i, l, s0, nmax, tol);
            return;
        else
            s = s + delta_s;
        end
    end

    fprintf('Convergence failed within %d iterations. Please try a different initial guess s0.\n', nmax);
end

function plot_shooting_results(i, t, sol, exact, Iter)
    shot_labels = "Shot " + string(1:i-1);
    figure;
    plot(t, sol, 'b', 'LineWidth', 2);
    hold on;
    plot(t, exact, 'r', 'LineWidth', 2);
    if i > 1
        plot(t, Iter', '--', 'LineWidth', 1);
        legend(['Numerical Solution', 'Exact Solution', shot_labels], 'Location', 'best');
    else
        legend({'Numerical Solution', 'Exact Solution'}, 'Location', 'best');
    end
    hold off;
    xlabel('t');
    ylabel('y(t)');
    title('Nonlinear BVP Solution via Newton-Raphson Shooting Method');
    grid on;
end

function export_shooting_table(t, sol, exact, error, i, l, s0, nmax, tol)
    u = length(t);
    RowNames = {'Exact', 'Approx'};

    fid = fopen('shooting_results.txt', 'a');
    fprintf(fid, '\n');
    fprintf(fid, 'Refinement level k:       %d\n', l);
    fprintf(fid, 'Initial slope guess s0:   %f\n', s0);
    fprintf(fid, 'Max allowed iterations:   %d\n', nmax);
    fprintf(fid, 'Tolerance:                %e\n', tol);
    fprintf(fid, 'Converged in shots:       %d\n', i);

    fprintf(fid, '+------------');
    for m = 1:u
        fprintf(fid, '+------------');
    end
    fprintf(fid, '+\n');

    fprintf(fid, '| %-10s ', 't');
    for m = 1:u
        fprintf(fid, '| %-10.5f ', t(m));
    end
    fprintf(fid, '|\n');

    fprintf(fid, '+------------');
    for m = 1:u
        fprintf(fid, '+------------');
    end
    fprintf(fid, '+\n');

    Data = [exact; sol];
    for m = 1:2
        fprintf(fid, '| %-10s ', RowNames{m});
        for p = 1:u
            fprintf(fid, '| %-10.5f ', Data(m, p));
        end
        fprintf(fid, '|\n');

        fprintf(fid, '+------------');
        for p = 1:u
            fprintf(fid, '+------------');
        end
        fprintf(fid, '+\n');
    end

    fprintf(fid, 'Maximum Absolute Error: %e\n', error);
    fclose(fid);
end

function f = ode_rhs(t, y)
    f = [y(2); -2 * y(1) * y(2) / t];
end

function g = variational_rhs(t, z, y)
    g = [z(2); (-2 / t) * (y(1) * z(2) + y(2) * z(1))];
end
