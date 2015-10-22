function [x] = newton ( x0, funder );
%--------------------------------------------------------------------------
%
% ... metodo de Newton rudimentario para ejemplificar el uso de la busqueda 
% lineal de M Overton. 
%
%       j-l morales ITAM
%       marzo 2015
%
%--------------------------------------------------------------------------
% ... tolerancias y valores iniciales
%
n = length(x0);
x = zeros(n,1);  F = zeros(n,1);  
%
TOL = 1.0d-8; k = 0; x = x0;  k = 0; ITMAX = 200;
c1 = 1.0e-4; c2 = 0.9; prtlevel = 2;

[ f, g, H ] = feval( funder, x );
%
norm_g = norm(g);

fprintf(1,'   k        f            ||g||    alpha     nfg \n');
fprintf(1,'------------------------------------------------\n');

while  norm_g > TOL  &  k < ITMAX
    p_N = -H\g;
    [ alpha, x, f, g, fail, nsteps ] = ...
                       linesch_sw(x, f, g, p_N, funder, c1, c2, prtlevel);
    
    [ ff, gg, H ] = feval( funder, x );
   
    norm_g   = norm(g);
    k = k + 1;
    fprintf(' %3i  %14.8e   %8.2e   %5.3f    %3i  \n', ...
                k,   f,      norm_g,  alpha, nsteps );   
end
