function stiff_ode_rlc_solver
    % Main program: Numerical integration of stiff and non-stiff linear ODE systems (RLC circuits)
    fprintf('--- Numerical Solver for Linear ODE Systems (RLC Circuits) ---\n');

    nprob = input('Select problem number (1: Non-stiff, 2: Stiff RLC system): ');
    nmet  = input('Select numerical method (1: Explicit Euler, 2: Explicit Trapezoidal, 3: Implicit Euler, 4: Implicit Trapezoidal): ');
    N     = input('Enter number of time steps N: ');

    [t0, tf, y0] = initial_conditions(nprob);
    tau = (tf - t0) / N;
    t = t0;
    y = y0;

    for n = 0:N-1
        ytau = step_method(nprob, nmet, t, y, tau);
        y = ytau;
        t = t + tau;
    end

    % Maximum absolute error at final time t = tf
    y_exact = exact_solution(nprob, tf);
    error_vec = abs(y - y_exact);
    max_error = max(error_vec);

    fprintf('\nMethod = %d | Problem = %d\n', nmet, nprob);
    fprintf('Steps N = %d | Max Error at final time (t = %.1f) = %13.4e\n', N, tf, max_error);
end

function [t0, tf, y0] = initial_conditions(nprob)
    t0 = 0;
    tf = 1;
    y0 = [0; 0];
    switch nprob
        case 1
            y0 = [4; 12];
        case 2
            y0 = [4; 18];
    end
end

function [A, b] = linear_system_matrix(nprob)
    switch nprob
        case 1
            R1 = 1;
            R2 = 3/5;
            If = 16;
            L  = 2;
            C  = 2/3;
        case 2
            R1 = 1;
            R2 = 2;
            If = 18;
            L  = 100;
            C  = 0.001;
    end
    A = [-R1/L, -1/L; 1/C, -1/(C*R2)];
    b = [R1/L * If; 0];
end

function f = derivative(nprob, ~, y)
    [A, b] = linear_system_matrix(nprob);
    f = A * y + b;
end

function y = exact_solution(nprob, t)
    switch nprob
        case 1
            lambda1 = -1;
            lambda2 = -2;
            a1 = -12; b1 = 6;  c1 = 10;
            a2 = -12; b2 = 18; c2 = 6;
        case 2
            lambda1 = -4.9997999879987992e+02;
            lambda2 = -3.0001200120011617e-02;
            a1 =  2.0002480392065366e-04;
            b1 = -2.0002000248039207e+00;
            c1 =  6;
            a2 =  1.0000640097615541e+01;
            b2 = -4.0006400976155421e+00;
            c2 =  12;
    end
    y = [a1 * exp(lambda1 * t) + b1 * exp(lambda2 * t) + c1;
         a2 * exp(lambda1 * t) + b2 * exp(lambda2 * t) + c2];
end

function ytau = step_method(nprob, nmet, t, y, tau)
    switch nmet
        case 1 % Explicit Euler
            ytau = y + tau * derivative(nprob, t, y);
        case 2 % Explicit Trapezoidal Rule (Heun's Method)
            f1 = derivative(nprob, t, y);
            f2 = derivative(nprob, t + tau, y + tau * f1);
            ytau = y + (tau / 2) * (f1 + f2);
        case 3 % Implicit Euler
            [A, b] = linear_system_matrix(nprob);
            Id = eye(size(A));
            if det(Id - tau * A) == 0
                disp('Error: Matrix (I - tau*A) is singular for the selected parameters.');
                ytau = y;
            else
                ytau = (Id - tau * A) \ (y + tau * b);
            end
        case 4 % Implicit Trapezoidal Rule (Crank-Nicolson for ODEs)
            [A, b] = linear_system_matrix(nprob);
            Id = eye(size(A));
            if det(Id - (tau / 2) * A) == 0
                disp('Error: Matrix (I - (tau/2)*A) is singular for the selected parameters.');
                ytau = y;
            else
                ytau = (Id - (tau / 2) * A) \ ((Id + (tau / 2) * A) * y + tau * b);
            end
    end
end
