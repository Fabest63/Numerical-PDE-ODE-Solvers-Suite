function fem_bspline_cubic
    % Finite Element Method (FEM) Solver using C^2 Cubic B-Splines (Order 3)
    mesh_list = [5, 10, 20, 40];

    a = 1;
    b = 2;
    x_eval = linspace(a, b, 101);
    u_exact = exact_solution(x_eval);

    for idx = 1:length(mesh_list)
        n = mesh_list(idx);
        h = (b - a) / n;
        xs = linspace(a, b, n + 1);

        % Assemble 7-band heptadiagonal sparse stiffness/mass matrix M
        di = zeros(n + 1, 1);
        for i = 1:n+1
            di(i) = integral(@(x) bilinear_form_integrand(x, i, i, n, h, a), ...
                max(a, xs(i) - 2*h), min(b, xs(i) + 2*h));
        end

        up1 = zeros(n + 1, 1);
        for i = 1:n
            up1(i + 1) = integral(@(x) bilinear_form_integrand(x, i, i + 1, n, h, a), ...
                max(a, xs(i) - h), min(b, xs(i) + 2*h));
        end
        lo1 = [up1(2:n+1); 0];

        up2 = zeros(n + 1, 1);
        for i = 1:n-1
            up2(i + 2) = integral(@(x) bilinear_form_integrand(x, i, i + 2, n, h, a), ...
                max(a, xs(i)), min(b, xs(i) + 2*h));
        end
        lo2 = [up2(3:n+1); 0; 0];

        up3 = zeros(n + 1, 1);
        for i = 1:n-2
            up3(i + 3) = integral(@(x) bilinear_form_integrand(x, i, i + 3, n, h, a), ...
                max(a, xs(i) + h), min(b, xs(i) + 2*h));
        end
        lo3 = [up3(4:n+1); 0; 0; 0];

        M = spdiags([lo3, lo2, lo1, di, up1, up2, up3], -3:3, n + 1, n + 1);

        % Assemble load vector v
        v = zeros(n + 1, 1);
        for i = 1:n+1
            v(i) = integral(@(x) load_vector_integrand(x, i, n, h, a), ...
                max(a, xs(i) - 2*h), min(b, xs(i) + 2*h));
        end

        % Solve sparse linear system for Cubic B-spline weights
        d = M \ v;

        % Evaluate Cubic B-Spline FEM solution on evaluation grid
        B = evaluate_bspline3_basis(n, h, a, x_eval);
        s = d' * B;

        max_error = max(abs(s - u_exact));

        figure;
        plot(x_eval, u_exact, '-r', 'LineWidth', 1.5, 'DisplayName', 'Exact Analytical Solution');
        hold on;
        plot(x_eval, s, '--b', 'LineWidth', 1.5, 'DisplayName', 'Cubic B-Spline FEM Approximation');
        grid on;
        legend('show', 'Location', 'best');
        title(sprintf('Sturm-Liouville FEM Solution (n = %d, Cubic B-Splines B_3)', n));
        xlabel('x'); ylabel('u(x)');
        hold off;

        fprintf('Maximum error with n = %2d (Cubic B-Spline):  %.4e\n', n, max_error);

        fid = fopen('fem_cubic_results.txt', 'a');
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

function b = master_bspline3(t)
    b = zeros(size(t));
    r1 = abs(t) < 1;
    r2 = abs(t) >= 1 & abs(t) <= 2;
    t1 = t(r1);
    t2 = t(r2);

    b(r1) = 2/3 - abs(t1).^2 + (abs(t1).^3) / 2;
    b(r2) = ((2 - abs(t2)).^3) / 6;
end

function B = evaluate_bspline3_basis(n, h, a, x)
    B = zeros(n + 1, length(x));
    t = (x - a) / h;
    Bm = master_bspline3(t + 1);
    Bp = master_bspline3(t - n - 1);
    for i = 1:n+1
        B(i, :) = master_bspline3(t - i + 1);
    end
    % Modify boundary splines to enforce homogeneous Dirichlet BCs
    B(2, :)   =  B(1, :) - 4 * B(2, :);
    B(1, :)   = -4 * Bm  + B(1, :);
    B(n, :)   = -4 * B(n, :) + B(n + 1, :);
    B(n+1, :) =  B(n + 1, :) - 4 * Bp;
end

function b = master_dbspline3(t)
    b = zeros(size(t));
    r1 = -2 <= t & t <= -1;
    r2 = -1 <  t & t <   0;
    r3 =  0 <= t & t <   1;
    r4 =  1 <= t & t <=  2;

    b(r1) =  0.5 * (2 + t(r1)).^2;
    b(r2) = -2 * t(r2) - 1.5 * t(r2).^2;
    b(r3) = -2 * t(r3) + 1.5 * t(r3).^2;
    b(r4) = -0.5 * (2 - t(r4)).^2;
end

function B = evaluate_dbspline3_basis(n, h, a, x)
    B = zeros(n + 1, length(x));
    t = (x - a) / h;
    Bm = (1 / h) * master_dbspline3(t + 1);
    Bp = (1 / h) * master_dbspline3(t - n - 1);
    for i = 1:n+1
        B(i, :) = (1 / h) * master_dbspline3(t - i + 1);
    end
    B(2, :)   =  B(1, :) - 4 * B(2, :);
    B(1, :)   = -4 * Bm  + B(1, :);
    B(n, :)   = -4 * B(n, :) + B(n + 1, :);
    B(n+1, :) =  B(n + 1, :) - 4 * Bp;
end

function y = bilinear_form_integrand(x, i, j, n, h, a)
    B  = evaluate_bspline3_basis(n, h, a, x);
    dB = evaluate_dbspline3_basis(n, h, a, x);
    y  = coeff_p(x) .* dB(i, :) .* dB(j, :) + coeff_q(x) .* B(i, :) .* B(j, :);
end

function y = load_vector_integrand(x, i, n, h, a)
    B = evaluate_bspline3_basis(n, h, a, x);
    y = source_f(x) .* B(i, :);
end



