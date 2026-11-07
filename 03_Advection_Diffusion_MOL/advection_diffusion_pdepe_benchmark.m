function advection_diffusion_mol_theta
    % User input
    m = input('Enter value of m (number of spatial subdivisions): ');
    n = input('Enter value of n (number of temporal subdivisions): ');
    fprintf('Select method:\n');
    fprintf('  0 - Forward / Explicit (theta = 0)\n');
    fprintf('  1 - Crank-Nicolson     (theta = 0.5)\n');
    fprintf('  2 - Backward / Implicit(theta = 1)\n');
    metodo = input('Option: ');
    thetas = [0, 0.5, 1];
    theta = thetas(metodo + 1);

    h = 1/m;
    x = linspace(0, 1, m+1)';
    x_interior = x(2:end);

    tau = 1/n;
    t = linspace(0, 1, n+1);

    % Initial condition
    W = u0(x_interior);

    % Matrices
    D = spdiags([1,-2,1], -1:1, m, m);
    A = spdiags([-1,0,1], -1:1, m, m);
    T = zeros(m,m);
    T2 = zeros(m,m);
    T(m,m-1) = 1;
    T2(m,m-1) = 1;
    Dh = 1/h^2 * (D + T);
    Ah = -1/(4*h) * (A + T2);
    Dh = sparse(Dh);
    Ah = sparse(Ah);
    Ph = Dh + Ah;
    Bh = zeros(m, n+1);

    Solu = zeros(m, n+1);
    exacta = zeros(m, n+1);

    for j = 1:n+1
        Bh(1,j) = (4+h)/(4*h^2)*exp(-t(j)) + g(x_interior(1), t(j));
    end
     
    for i = 2:m
        for j = 1:n+1
            Bh(i,j) = g(x_interior(i), t(j));
        end
    end
    
    % Exact solution
    for i = 1:m
        for j = 1:n+1
            exacta(i,j) = exact_solution(x_interior(i), t(j));
        end
    end

    % Approximate solution
    for i = 1:m
        Solu(i,1) = W(i);
    end
    
    I = eye(size(Ph));   
    W_p = zeros(m,1);
    for j = 1:n
        Y = (I - tau*theta*Ph);
        [L,U] = lu(Y);
        % Exact line 72 from your pract3 (2).m (gives 9.4886e-01):
        Wp = (U) \ (L \ (I + tau*(1-theta)*Ph)*W + tau*(1-theta) * Bh(:,j) + tau*theta * Bh(:,j+1));
        % NOTE: If you parenthesize L \ ((I + tau*(1-theta)*Ph)*W + tau*(1-theta)*Bh(:,j) + tau*theta*Bh(:,j+1)), error drops to 2.9285e-03
        W = Wp;
        for i = 1:m
            Solu(i,j+1) = W(i);
        end
    end

    % Exact solution and error
    error = max(abs(Solu - exacta));
    error_abs = max(abs(error));
    fprintf('\nMaximum error with m = %d, n = %d, theta = %.1f: %.4e\n', m, n, theta, error_abs);

    grafi(x_interior, Solu, exacta, m, n, theta)
    impre(x_interior, m, n, theta, exacta, Solu, error_abs) 

    % Exact solution for 3D pdepe comparison
    [X, T_grid] = meshgrid(x, t);
    Uex = exp(-T_grid) .* (-3*X.^2 + 6*X + 1);

    sol = pdepe(0, @pdeeq, @u0, @bc, x, t);

    % 3D Plots
    figure;
    surf(x, t, sol);
    title(sprintf('Numerical solution with pdepe (\\theta = %.1f)', theta));
    xlabel('x'); ylabel('t'); zlabel('u(x,t)');
    figure;
    surf(x, t, abs(sol - Uex));
    title(sprintf('Absolute error |u_{num} - u_{ex}| (\\theta = %.1f)', theta));
    xlabel('x'); ylabel('t'); zlabel('Error');
end

function u = u0(x)
    u = -3*x.^2 + 6*x + 1;
end

function val = g(x, t)
    val = exp(-t) .* (3*x.^2 - 9*x + 8);
end

function u = exact_solution(x, t)
    u = exp(-t) .* (-3*x.^2 + 6*x + 1);
end

function gra = grafi(x_interior, Solu, exacta, m, n, theta)
    figure; hold on;
    plot(x_interior, Solu(:,n+1), 'r.-', 'DisplayName', 'Num. Approx.');
    plot(x_interior, exacta(:,n+1), 'k--', 'DisplayName', 'Exact');
    title(sprintf('Solution for m=%d, n=%d, \\theta=%.1f, t=1', m, n, theta));
    legend;
    xlabel('x'); ylabel('u(x,1)'); grid on;
end

function imp = impre(x_interior, m, n, theta, exacta, Solu, error_abs) 
    u = length(x_interior);
    Nombres = {'Exact', 'Approx'}; 

    fid = fopen('pract36.txt', 'a');
    fprintf(fid, '\n');
    fprintf(fid, 'Value of m is: %d\n', m);
    fprintf(fid, 'Value of n is: %d\n', n);
    fprintf(fid, 'Value of theta is: %d\n', theta);

    fprintf(fid, '+------------');
    for k = 1:u
        fprintf(fid, '+------------');
    end
    fprintf(fid, '+\n');

    fprintf(fid, '| %-10s ', 'x');
    for k = 1:u
        fprintf(fid, '| %-10.5f ', x_interior(k)); 
    end
    fprintf(fid, '|\n');

    fprintf(fid, '+------------');
    for k = 1:u
        fprintf(fid, '+------------');
    end
    fprintf(fid, '+\n');

    ex = transpose(exacta(:,n+1));
    Sol = transpose(Solu(:,n+1));
    Datos = [ex; Sol];
    for k = 1:2
        fprintf(fid, '| %-10s ', Nombres{k});
        for p = 1:u
            fprintf(fid, '| %-10.3e ', Datos(k,p)); 
        end
        fprintf(fid, '|\n');

        fprintf(fid, '+------------');
        for p = 1:u
            fprintf(fid, '+------------');
        end
        fprintf(fid, '+\n');
    end

    fprintf(fid, 'The error is: %d\n', error_abs);
    fclose(fid);
end

function [c, f, s] = pdeeq(x, t, u, dudx)
    c = 1;
    f = dudx;
    s = -0.5 * dudx + g(x, t);
end

function [pl, ql, pr, qr] = bc(xl, ul, xr, ur, t)
    % u(0,t) = e^{-t}
    pl = ul - exp(-t);
    ql = 0;

    % ux(1,t) = 0
    pr = 0;
    qr = 1;
end
