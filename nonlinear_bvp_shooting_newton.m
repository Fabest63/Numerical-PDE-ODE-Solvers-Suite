function nonlinear_bvp_shooting_newton
l = input('Ingrese un valor de k entre 2 y 6: ');
s0 = input('Ingrese un valor para s0: ');
nmax = input('Ingrese el número máximo de iteraciones: ');
tol = input('Ingrese la tolerancia: ');
valores_l = [2, 3, 4, 5, 6];
if ~ismember(l, valores_l) 
    disp('No es válido');
else
    alpha = 1/2;
    beta = 2/3;
    a = 1;
    b = 2;
    n = 2^l;
    tau = 1/n;
    t = linspace(a,b,n+1);
    s = s0;
    i = 0;
    Iter = zeros(nmax,n+1);
    while i <= nmax-1
        i=i+1;
        % Resolucion del PVI(s)
        y0 = [alpha;s];
        y = y0;
        Y = zeros(2,n+1);
        Y(:,1) = y0;
        for k=0:n-1
            f1 = der(t(k+1),y);
            y1 = y+tau*f1;
            f2 =der(t(k+2),y1);
            ytau = y + tau*(f1+f2)/2;
            Y(:,k+2)=ytau;
            y = ytau; 
        end

        %Asi w(b,s) = yn,1 (columna última de Y fila 1)

        %Resolucion del PVI(xi)
        z0 = [0;1];
        z = z0;
        Z = zeros(2,n+1); %No es necesario
        Z(:,1) = z0; %No es necesario
        for j=0:n-1
            g1 = der2(t(j+1),z,Y(:,j+1));
            z1 = z+tau*g1;
            g2 = der2(t(j+2),z1,Y(:,j+2));
            ztau = z + tau*(g1+g2)/2;
            Z(:,j+2)=ztau; %No es necesario
            z = ztau;
        end

        difer = -(Y(1,n+1)-beta)/z(1);
        Iter(i,:) = Y(1,:);
        Iter(all(abs(Iter) < tol, 2), :) = [];
        

        %Control de parada

        if abs(difer) <= tol
            Iter(end, :) = []; 
            sol = Y(1,:); %Aproximacion a w(t_k)
            exacta = t./(1+t);
            error = max(abs(sol-exacta));
            disp(['El valor de la solución exacta en los puntos de aproximación es ', num2str(exacta)]);
            disp(['El valor de la solución aproximada en los puntos de aproximación es: ', num2str(sol)]);
            disp(['El error es: ', num2str(error)]);
            disp(['El número de tiros es de: ', num2str(i)]);

           grafi(i,t,sol,exacta,Iter)
           impre(t,sol,exacta,error,i,l,s0,nmax,tol)

            
            return
        else
            s = s + difer;

        end

    end

    fprintf("Fallo de convergenica, introduzca otro s0")

end
end


function graf = grafi(i,t,sol,exacta,Iter)
num_tiros = "Tiro " + string(1:i-1);
figure;
plot(t, sol, 'b', 'LineWidth', 2); 
hold on; 
plot(t, exacta, 'r', 'LineWidth', 2); 
hold on;
plot(t, Iter', '--','LineWidth', 1); 
hold off; 
xlabel('t');
ylabel('Valor de la función');
title('Gráfico de la solución exacta y la aproximada');
legend('Aproximada', 'Exacta');
legend(['Solución aproximada','Solución exacta', num_tiros, ], 'Location', 'best');
grid on;
end


function imp = impre(t,sol,exacta,error,i,l,s0,nmax,tol)
           u = length(t); % Número de columnas

            % Definir nombres de las filas
            Nombres = {'Exacta', 'Aproximada'}; 

            % Crear archivo de texto
            fid = fopen('shooting_results.txt', 'a');
            fprintf(fid, '\n');

            fprintf(fid, 'El valor de k es: %d\n', l);
            fprintf(fid, 'El valor de s0 es: %d\n', s0);
            fprintf(fid, 'El número máximo de iteraciones es de: %d\n', nmax);
            fprintf(fid, 'La tolerancia es de: %d\n', tol);
            fprintf(fid, 'El número de tiros fue de: %d\n', i);

            % Imprimir línea superior de la tabla
            fprintf(fid, '+------------');
            for m = 1:u
               fprintf(fid, '+------------');
            end
            fprintf(fid, '+\n');

            % Imprimir encabezados (Primera celda vacía + valores del Tiempo)
            fprintf(fid, '| %-10s ', 't'); % Espacio para la primera columna
            for m = 1:u
                fprintf(fid, '| %-10.5f ', t(m)); 
            end
            fprintf(fid, '|\n');

            % Imprimir línea de separación
            fprintf(fid, '+------------');
            for m = 1:u
               fprintf(fid, '+------------');
            end
            fprintf(fid, '+\n');

            % Imprimir datos con alineación centrada
            Datos = [exacta; sol]; % Matriz con los valores
            for m = 1:2
                fprintf(fid, '| %-10s ', Nombres{m}); % Nombre de fila
                for p = 1:u
                    fprintf(fid, '| %-10.5f ', Datos(m, p)); 
                end
                fprintf(fid, '|\n');

                % Imprimir línea de separación después de cada fila
                fprintf(fid, '+------------');
                for p = 1:u
                   fprintf(fid, '+------------');
                end
                fprintf(fid, '+\n');
            end

        fprintf(fid, 'El error es: %d\n', error);

        fclose(fid);

end


function f=der(t,y)
f=[y(2);-2*y(1)*y(2)/t];
end


function g=der2(t,z,y)
g=[z(2);-2/t*(y(1)*z(2)+y(2)*z(1))];
end

