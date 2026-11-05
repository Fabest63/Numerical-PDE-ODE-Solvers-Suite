% Metodo de lineas para la ecuacion de adveccion-difusion

function advection_diffusion_mol_theta

    % Entrada del usuario
    m = input('Ingrese el valor de m (nro de divisiones espaciales): ');
    n = input('Ingrese el valor de n (nro de divisiones temporales): ');
    fprintf('Seleccione el metodo:\n');
    fprintf('  0 - Explicito (theta = 0)\n');
    fprintf('  1 - Crank-Nicolson (theta = 0.5)\n');
    fprintf('  2 - Implicito (theta = 1)\n');
    metodo = input('Opcion: ');
    thetas = [0, 0.5, 1];
    theta = thetas(metodo + 1);

    h = 1/m;
    x = linspace(0, 1, m+1)';
    x_interior = x(2:end);

    tau = 1/n;
    t = linspace(0, 1, n+1);

    % Condicion inicial
    W = u0(x_interior);
    W0 = W;

    % Matrices
    [D, B] = construir_matrices(m, h);
    D
    I = speye(m);

    % Paso en el tiempo
    for j = 1:n
        tj = t(j);
        tj1 = t(j+1);

        bhj = g(x_interior, tj);
        bhj1 = g(x_interior, tj1);

        rhs = (I + tau*(1-theta)*D)*W + tau*(1-theta)*bhj + tau*theta*bhj1;
        A = (I - tau*theta*D);
        [L, U] = lu(A);
        W = U \ (L \ rhs);
    end

    % Solución exacta y error
    u_ex = exact_solution(x_interior, 1);
    error = max(abs(W - u_ex));
    fprintf('\nError maximo con m = %d, n = %d, theta = %.1f: %.4e\n', m, n, theta, error);

    % Graficar
    figure; hold on;
    plot(x_interior, W, 'r.-', 'DisplayName', 'Aprox. num.');
    plot(x_interior, exact_solution(x_interior, 1), 'k--', 'DisplayName', 'Exacta');
    title(sprintf('Solución para m=%d, n=%d, \theta=%.1f, t=1', m, n, theta));
    legend;
    xlabel('x'); ylabel('u(x,1)'); grid on;
end

% ------------------------------
% FUNCIONES AUXILIARES
% ------------------------------

function u = u0(x)
    u = -3*x.^2 + 6*x + 1;
end

function val = g(x, t)
    val = exp(-t) .* (3*x.^2 - 9*x + 8);
end

function u = exact_solution(x, t)
    u = exp(-t) .* (-3*x.^2 + 6*x + 1);
end

function [D, B] = construir_matrices(m, h)
    % D es la matriz para el operador uxx - (1/2)ux
    e = ones(m, 1);
    % Derivada segunda central
    Axx = spdiags([e -2*e e], -1:1, m, m) / h^2;
    % Derivada primera central
    Ax = spdiags([-e zeros(m,1) e], -1:1, m, m) / (2*h);
    D = Axx - 0.5 * Ax;

    % Condiciones de frontera
    B = zeros(m,1);
    B(1) = B(1) + (1/h^2 + 0.5/(2*h)) * exp(-0); % u(0,t) = e^{-t}, t=0
end
