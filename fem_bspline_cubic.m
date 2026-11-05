function fem_bspline_cubic

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
di = zeros(n+1,1);
for i = 1:n+1
    di(i) = integral(@(x)masa(x,i,i,n,h,a),...
        max(a,xs(i) - 2*h),min(b,xs(i) + 2*h));
end

up1 =zeros(n+1,1);
for i=1:n
    up1(i+1) = integral(@(x)masa(x,i,i+1,n,h,a),...
        max(a,xs(i) - h),min(b,xs(i) + 2*h));
end
lo1 = [up1(2:n+1);0];

up2 =zeros(n+1,1);
for i=1:n-1
    up2(i+2) = integral(@(x)masa(x,i,i+2,n,h,a),...
        max(a,xs(i)),min(b,xs(i) + 2*h));
end
lo2 = [up2(3:n+1);0;0];

up3 =zeros(n+1,1);
for i=1:n-2
    up3(i+3) = integral(@(x)masa(x,i,i+3,n,h,a),...
        max(a,xs(i) + h),min(b,xs(i) + 2*h));
end
lo3 = [up3(4:n+1);0;0;0];

M = spdiags([lo3,lo2,lo1,di,up1,up2,up3],-3:3,n+1,n+1);
% Vector v
v = zeros(n+1,1);
for i = 1:n+1
    v(i) = integral(@(x)rhs(x,i,n,h,a),...
        max(a,xs(i)-2*h),min(b,xs(i) + 2*h));
end

d = M\v;

B = bspline3(n,h,a,x_error);

s = d'*B;

error = max(abs(s-u));

figure;
plot(x_error,u,"-r","DisplayName","Solución real")
hold on
plot(x_error,s,"-b","DisplayName","Solución Aproximada")
grid on
legend("show");
title(sprintf('Solución para n=%d con spline de orden 3', n));
xlabel("Eje X")
hold off

fprintf('EL error maximo con n=%d es: %.4e\n', n,  error);

     % Crear archivo de texto
     fid = fopen('fem_cubic_results.txt', 'a');
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


function b=patronB3(t)

b = zeros(size(t));

rango1 = abs(t)< 1;
rango2 = abs(t) >= 1 & abs(t) <=2;
t1 = t(rango1);
t2 = t(rango2);

b(rango1) = 2/3 - abs(t1).^2 + abs(t1) .^3 / 2;
b(rango2) = (2-abs(t2)).^3/6;

end


function B = bspline3(n,h,a,x)

B = zeros(n+1,length(x));
t = (x-a)/h;
Bm = patronB3(t+1);
Bp = patronB3(t-n-1);
for i=1:n+1
    B(i,:) = patronB3(t-i+1);
end
B(2,:) = B(1,:) - 4*B(2,:); 
B(1,:) = -4*Bm + B(1,:);
B(n,:) = -4*B(n,:) + B(n+1,:);
B(n+1,:) = B(n+1,:) - 4*Bp;
end

function b=patrondB3(t)

b = zeros(size(t));

rango1 = -2 <= t & t <= -1;
rango2 = -1 < t & t < 0;
rango3 = 0 <= t & t < 1;
rango4 = 1 <= t & t <= 2;
t1 = t(rango1);
t2 = t(rango2);
t3 = t(rango3);
t4 = t(rango4);

b(rango1) = 1/2 .* (2+t1).^2;
b(rango2) = -2.*t2 -3/2 .* t2.^2;
b(rango3) = -2.*t3 + 3/2 .* t3.^2;
b(rango4) = -1/2 .* (2-t4).^2;

end

function B = dbspline3(n,h,a,x)

B = zeros(n+1,length(x));
t = (x-a)/h;
Bm = 1/h * patrondB3(t+1);
Bp = 1/h * patrondB3(t-n-1);
for i=1:n+1
    B(i,:) = 1/h *patrondB3(t-i+1);
end
B(2,:) = B(1,:) - 4*B(2,:); 
B(1,:) = -4*Bm + B(1,:);
B(n,:) = -4*B(n,:) + B(n+1,:);
B(n+1,:) = B(n+1,:) - 4*Bp;
end



function y = masa(x,i,j,n,h,a)

B = bspline3(n,h,a,x);
dB = dbspline3(n,h,a,x); 
y = p(x) .*dB(i,:) .*dB(j,:) + q(x) .*B(i,:) .*B(j,:);

end

function y = rhs(x,i,n,h,a)

  B = bspline3(n,h,a,x); 
  y = f(x) .*B(i,:);
  
end



