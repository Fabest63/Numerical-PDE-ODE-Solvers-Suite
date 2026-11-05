function fem_bspline_linear

lista = [5,10,20,40];
l = 1;

while l < length(lista) + 1
n = lista(l);
l = l + 1;
a = 1;
b = 2;
h = (b-a)/n;
xs = linspace(a,b,n+1);

%Exacta
x_error=linspace(a,b,101);
u = exacta(x_error);

%Matriz masa
di = zeros(n-1,1);
for i = 1:n-1
    di(i) = integral(@(x)masa(x,i,i,n,h,a),...
        xs(i),xs(i) + 2*h);
end

up1 =zeros(n-1,1);
for i=1:n-2
    up1(i+1) = integral(@(x)masa(x,i,i+1,n,h,a),...
        xs(i) + h,xs(i) + 2*h);
end
lo1 = [up1(2:n-1);0];


M = spdiags([lo1,di,up1],-1:1,n-1,n-1);
% Vector v
v = zeros(n-1,1);
for i = 1:n-1
    v(i) = integral(@(x)rhs(x,i,n,h,a),...
        xs(i),xs(i) + 2*h);
end

d = M\v;

B = bspline1(n,h,a,x_error);

s = d'*B;

error = max(abs(s-u));

figure;
plot(x_error,u,"-r","DisplayName","Solución real")
hold on
plot(x_error,s,"-b","DisplayName","Solución Aproximada")
grid on
legend("show");
title(sprintf('Solución para n=%d con spline de orden 1', n));
xlabel("Eje X")
hold off

fprintf('EL error maximo con n=%d es: %.4e\n', n,  error);

     % Crear archivo de texto
     fid = fopen('fem_linear_results.txt', 'a');
     fprintf(fid, '\n');

     fprintf(fid, 'El valor de n es: %d\n', n);
     fprintf(fid, 'El error es: %d\n', error);
     fclose(fid);
end
end


function y=p(x)
y = x.^2;
end
function y=q(x) 
y = 6;
end

function y=f(x) 
y = 6 * x.^3;
end

function u = exacta(x)
u = 63/31 * x.^2 - 32./(31*x.^3) - x.^3;
end

function b=patronB1(t)

b = zeros(size(t));

rango1 = abs(t)< 1;
t1 = t(rango1);

b(rango1) = 1 - abs(t1);


end


function B = bspline1(n,h,a,x)

B = zeros(n-1,length(x));
t = (x-a)/h;
for i=1:n-1
    B(i,:) = patronB1(t-i);
end
end

function b=patrondB1(t)

b = zeros(size(t));

rango1 = -1 < t & t < 0;
rango2 = 0 < t & t < 1;


b(rango1) = 1;
b(rango2) = -1;


end

function B = dbspline1(n,h,a,x)

B = zeros(n-1,length(x));
t = (x-a)/h;

for i=1:n-1
    B(i,:) = 1/h *patrondB1(t-i);
end

end


function y = masa(x,i,j,n,h,a)

    B = bspline1(n,h,a,x);
    dB = dbspline1(n,h,a,x); 
    y = p(x) .*dB(i,:) .*dB(j,:) + q(x) .*B(i,:) .*B(j,:);

end

function y = rhs(x,i,n,h,a)

 B = bspline1(n,h,a,x);
  y = f(x) .*B(i,:);

end



