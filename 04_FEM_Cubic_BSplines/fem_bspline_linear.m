function fem_bspline_linear
    % Finite Element Method (FEM) Solver using Linear B-Splines (Order 1)
    mesh_list = [5, 10, 20, 40];

    a = 1;
    b = 2;
    x_eval = linspace(a, b, 101);
    u_exact = exact_solution(x_eval);

    for idx = 1:length(mesh_list)
        n = mesh_list(idx);
        h = (b - a) / n;
        xs = linspace(a, b, n + 1);

        % Assemble tridiagonal stiffness/mass matrix M via adaptive quadrature
        di = zeros(n - 1, 1);
        for i = 1:n-1
            di(i) = integral(@(x) bilinear_form_integrand(x, i, i, n, h, a), xs(i), xs(i) + 2 * h);
        end

        up1 = zeros(n - 1, 1);
        for i = 1:n-2
            up1(i + 1) = integral(@(x) bilinear_form_integrand(x, i, i + 1, n, h, a), xs(i) + h, xs(i) + 2 * h);
        end
        lo1 = [up1(2:n-1); 0];

        M = spdiags([lo1, di, up1], -1:1, n - 1, n - 1);

        % Assemble load vector v
        v = zeros(n - 1, 1);
        for i = 1:n-1
            v(i) = integral(@(x) load_vector_integrand(x, i, n, h, a), xs(i), xs(i) + 2 * h);
        end

        % Solve sparse linear system for B-spline coefficients
        d = M \ v;

        % Evaluate FEM approximation on fine grid
        B = evaluate_bspline1_basis(n, h, a, x_eval);
        s = d' * B;

        max_error = max(abs(s - u_exact));

        figure;
        plot(x_eval, u_exact, '-r', 'LineWidth', 1.5, 'DisplayName', 'Exact Analytical Solution');
        hold on;
        plot(x_eval, s, '--b', 'LineWidth', 1.5, 'DisplayName', 'Linear B-Spline FEM Approximation');
        grid on;
        legend('show', 'Location', 'best');
        title(sprintf('Sturm-Liouville FEM Solution (n = %d, Linear B-Splines B_1)', n));
        xlabel('x'); ylabel('u(x)');
        hold off;

        fprintf('Maximum error with n = %2d (Linear B-Spline): %.4e\n', n, max_error);

        fid = fopen('fem_linear_results.txt', 'a');
        fprintf(fid, '\nMesh subdivisions n = %d | Max Error = %.6e\n', n, max_error);
        fclose(fid);
    end
end

function y = coeff_p(x); y = x.^2; end
function y = coeff_q(~); y = 6; end
function y = source_f(x); y = 6 * x.^3; end

function u = exact_solution(x)
    u = (63 / 31) * x.^2 - 32 ./ (31 * x.^3) - x.^3;
end

function b = master_bspline1(t)
    b = zeros(size(t));
    mask = abs(t) < 1;
    b(mask) = 1 - abs(t(mask));
end

function B = evaluate_bspline1_basis(n, h, a, x)
    B = zeros(n - 1, length(x));
    t = (x - a) / h;
    for i = 1:n-1
        B(i, :) = master_bspline1(t - i);
    end
end

function b = master_dbspline1(t)
    b = zeros(size(t));
    b(-1 < t & t < 0) =  1;
    b( 0 < t & t < 1) = -1;
end

function B = evaluate_dbspline1_basis(n, h, a, x)
    B = zeros(n - 1, length(x));
    t = (x - a) / h;
    for i = 1:n-1
        B(i, :) = (1 / h) * master_dbspline1(t - i);
    end
end

function y = bilinear_form_integrand(x, i, j, n, h, a)
    B  = evaluate_bspline1_basis(n, h, a, x);
    dB = evaluate_dbspline1_basis(n, h, a, x);
    y  = coeff_p(x) .* dB(i, :) .* dB(j, :) + coeff_q(x) .* B(i, :) .* B(j, :);
end

function y = load_vector_integrand(x, i, n, h, a)
    B = evaluate_bspline1_basis(n, h, a, x);
    y = source_f(x) .* B(i, :);
end



