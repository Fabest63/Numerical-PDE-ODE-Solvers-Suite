% Resolución del problema de advección-difusión usando PDEPE

function advection_diffusion_pdepe_benchmark

    % Entrada del usuario
    m = input('Ingrese el valor de m (nro de divisiones espaciales): ');
    n = input('Ingrese el valor de n (nro de divisiones temporales): ');

    fprintf('Seleccione el metodo para comparar:\n');
    fprintf('  0 - Explicito (theta = 0)\n');
    fprintf('  1 - Crank-Nicolson (theta = 0.5)\n');
    fprintf('  2 - Implicito (theta = 1)\n');
    metodo = input('Opcion: ');
    thetas = [0, 0.5, 1];
    theta = thetas(metodo + 1);

    x = linspace(0, 1, m+1);
    t = linspace(0, 1, n+1);

    sol = pdepe(0, @pdeeq, @u0, @bc, x, t);

    % Solucion exacta
    [X, T] = meshgrid(x, t);
    Uex = exp(-T) .* (-3*X.^2 + 6*X + 1);

    % Error
    error = max(max(abs(sol - Uex)));
    fprintf('Error maximo con pdepe (comparado con \theta = %.1f): %.4e\n', theta, error);

    % Grafica 3D
    figure;
    surf(x, t, sol);
    title(sprintf('Solución numérica con pdepe (\theta = %.1f)', theta));
    xlabel('x'); ylabel('t'); zlabel('u(x,t)');

    figure;
    surf(x, t, abs(sol - Uex));
    title(sprintf('Error absoluto |u_{num} - u_{ex}| (\theta = %.1f)', theta));
    xlabel('x'); ylabel('t'); zlabel('Error');
end

% ------------------------------
% FUNCIONES AUXILIARES
% ------------------------------

function [c, f, s] = pdeeq(x, t, u, dudx)
    c = 1;
    f = dudx;
    s = -0.5 * dudx + g(x, t);
end

function u = u0(x)
    u = -3*x.^2 + 6*x + 1;
end

function [pl, ql, pr, qr] = bc(xl, ul, xr, ur, t)
    % u(0,t) = e^{-t}
    pl = ul - exp(-t);
    ql = 0;

    % ux(1,t) = 0
    pr = 0;
    qr = 1;
end

function val = g(x, t)
    val = exp(-t) .* (3*x.^2 - 9*x + 8);
end
