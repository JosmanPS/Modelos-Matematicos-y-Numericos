function [ x, normg, iter, t ] = BFGS_BLOverton(funAMPL, maxiter )
%Alejandra Lelo de Larrea Ibarra
%1244
%�Seguimos usando diferencias finitas? o los valores que arroja linesch_sw???

nombre_ampl = strcat(funAMPL, '.nl');

[ x0, xlow, xupp, y0, clow, cupp ] = spamfunc( nombre_ampl );
n = length(x0);
x = x0;

% [ f, c ] = spamfunc( x, 0 );                  % evaluar la funcion objetivo
% [ g, A ] = spamfunc( x, 1 );                  % evaluar el gradiente 
% [ H ]    = spamfunc( y0 ); 

tol=10^-6;
n=length(x0);

tic;

iter = 0;
H = eye(n);
[f, c] = spamfunc(x, 0);
[g, A] = spamfunc(x, 1);
normg = norm(g);
curv = 1;
while normg > tol && maxiter > iter && curv > 0
    
    p = -H * g;
    [alpha, xx, ff, gg, fail, nsteps] = linesch_swAMPL(x, f, g, p, funAMPL, 10^-4, 0.9, 0);
    s = xx - x;
    y=gg-g;
    curv = s'*y;
    H = H + ((s'*y+y'*H*y)*(s*s'))/((s'*y)^2)-((H*y*s'+s*y'*H)/(s'*y));
    normg = norm(gg,2)
    x = xx;
    g = gg;
    f = ff;
    iter = iter + 1;
    
end

t = toc;

end

