function [x] = BFGS ( x0, funder, maxiter, tol );
%
% metodo de Newton global (busqueda lineal con condiciones fuertes
% de Wolfe) para optimizacion de funciones sin restricciones
% j-l morales 19 de febrero, 2015, ITAM
%
%--------------------------------------------------------------------------
%
% ... tolerancias y valores iniciales
%
n = length(x0);
x = x0;  
k = 0; 
%
maxgc = 2*n;
tolgc = tol;
c1 = 1.0e-4; c2 = 0.9;

f = feval( funder, x );
g = diffgrad( x, funder );
H = eye(n);

norm_g = norm(g);

fprintf(1,'   k        f            ||g||     alfa     nfg     curv  \n');
fprintf(1,'-----------------------------------------------------------------------\n');

while  norm_g > tol  &&  k < maxiter
    
    x1 = x;
    g1 = g;
    
    p_N = - H * g;
    [ alfa, x, f, ~, ~, numfg ] = biseccion ( x, f, g, p_N, funder, c1, c2 );

    g = diffgrad( x, funder );
    
    s = x - x1;
    y = g - g1;
    
    H = H + (s'*y + y'*H*y)*(s*s')/(s'*y)^2 - (H*y*s' + s*y'*H)/(s'*y);
 
    curv = s' * y;
    
    norm_g   = norm(g); 
    k = k + 1;
    fprintf(' %3i  %14.8e   %8.2e   %5.3f    %3i       %5.3f   \n', ...
                            k, f,  norm_g, alfa, numfg, curv );
end
%
%==========================================================================
%