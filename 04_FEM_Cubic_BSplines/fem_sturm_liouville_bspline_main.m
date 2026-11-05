function fem_sturm_liouville_bspline_main

 
n = input("Introduce si quieres hacerlo por B-spline de orden 1 o de orden 3 (1 o 3): ");
if n == 1
    fem_bspline_linear
elseif n == 3
    fem_bspline_cubic
else
    fprintf("Error en la selección del método por favor, elija uno válido.")
end

end

