function stiff_ode_rlc_solver

fprintf("Este programa intera problemas\n");

nprob = input("Introduce el número de problema 1 o 2: ");
nmet = input("Introduce el número de metodo 1 a 4: ");
N = input("Introduce el número de pasos: ");
[t0,tf,y0] = inicial(nprob);
tau = (tf-t0)/N;
t = t0;
y = y0;
for n=0:N-1
    ytau = metodo(nprob,nmet,t,y,tau); 
    y = ytau;
    t = t + tau;
end

%Error en el punto final
yexacta = solucion(nprob,tf);
error = abs(y-yexacta);
maxerror = max(error);
fprintf("Metodo= %d,Problema= %d \n ", nmet,nprob);
fprintf("N= %d, Error en el punto final= %13.4e \n",N, maxerror);
end



function [t0,tf,y0] = inicial(nprob)
t0 = 0;
tf= 1;
y0 = [0;0];
switch nprob
    case 1
    y0 = [4;12];
    case 2
    y0 = [4;18];
end
end



function [A,b] = lin(nprob)
switch nprob

    case 1
        R1=1;
        R2=3/5;
        If=16;
        L=2;
        C=2/3;
    case 2
        R1=1;
        R2=2;
        If=18;
        L=100;
        C=0.001;
end
A=[-R1/L, -1/L; 1/C, -1/(C*R2)];
b=[R1/L * If;0];
end



function f=derivada(nprob,t,y)
[A,b]=lin(nprob);
f = A*y+b;
end



function y = solucion(nprob,t)
%syms t
%f = derivada(nprob,t,y);
%y = int(f);

switch nprob
    case 1
        lambda1 = -1;
        lambda2 = -2;
        a1 = -12;
        b1 = 6;
        c1 = 10;
        a2 = -12;
        b2 = 18;
        c2 = 6;
    case 2
        lambda1 = -4.9997999879987992e+02;
        lambda2 = -3.0001200120011617e-02;
        a1 = 2.0002480392065366e-04;
        b1 = -2.0002000248039207e+00;
        c1 = 6;
        a2 = 1.0000640097615541e+01;
        b2 = -4.0006400976155421e+00;
        c2 = 12;
end
y = [a1*exp(lambda1*t)+b1*exp(lambda2*t)+c1; 
    a2*exp(lambda1*t)+b2*exp(lambda2*t)+c2];
end



function ytau=metodo(nprob,nmet,t,y,tau)
switch nmet
    case 1 %Euler explicito
        ytau = y + tau+derivada(nprob,t,y);
    case 2 %Regla trapezoidal explicita
        f1 = derivada(nprob,t,y);
        f2 = derivada(nprob,t+tau,y+tau *f1);
        ytau = y + tau/2 * f1 + tau/2 * f2;
    case 3 %Euler implicito
        [A,b] = lin(nprob);
        Id =eye(size(A));
        if det(Id-tau*A) == 0
            disp("El método no funciona con los valores introducidos la matriz no es invertible")
        else
            ytau = (Id-tau*A)\(y+tau*b);
        end
end
end