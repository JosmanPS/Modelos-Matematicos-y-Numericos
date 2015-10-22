function [ alfa, x, f, g, falla, numfg ] = biseccion (x, f, g, p, funder, c1, c2 );
%
% busqueda lineal simplificada 
% j-l morales 19 defebrero, 2015, ITAM
%
%--------------------------------------------------------------------------

falla = 0;                 %  1 si el procedimiento hizo mas de 20 pasos
numfg = 0;                 %  numero de evaluaciones de f y g
alfa  = 1;  gTp = g'*p;    %  derivada direccional
%
while numfg < 20
    xp = x + alfa*p;
    [ fp, gp, H ] = feval( funder, xp );
    numfg = numfg + 1;
    if fp <= f + c1*alfa*gTp & abs(gp'*p) <= c2*abs(gTp)
       x = xp; f = fp; g = gp;
       return
    else
        alfa = alfa/2;
    end
end
falla = 1;
% 
% 
% newton([-1;1],@rosenbrock, 10000, 1.0e-6);
%    k        f            ||g||     alfa     nfg     GC           ||r|| 
% -----------------------------------------------------------------------
%    1  3.45312500e+00   2.55e+01   0.125      4  residuo < tol  0.00e+00 
%    2  2.65383704e+00   8.14e+00   1.000      1  residuo < tol  6.43e-14 
%    3  2.24223200e+00   1.34e+01   0.500      2  residuo < tol  2.58e-13 
%    4  1.70115714e+00   6.81e+00   1.000      1  residuo < tol  1.32e-14 
%    5  1.40547604e+00   1.14e+01   1.000      1  residuo < tol  4.99e-15 
%    6  9.33002749e-01   2.37e+00   1.000      1  residuo < tol  8.01e-15 
%    7  7.56689459e-01   8.25e+00   0.500      2  residuo < tol  4.89e-15 
%    8  4.77302567e-01   1.55e+00   1.000      1  residuo < tol  5.26e-15 
%    9  3.53089573e-01   5.49e+00   0.500      2  residuo < tol  2.76e-15 
%   10  2.11933410e-01   2.30e+00   1.000      1  residuo < tol  3.89e-14 
%   11  1.43174352e-01   7.62e+00   1.000      1  residuo < tol  1.73e-14 
%   12  6.06288353e-02   6.38e-01   1.000      1  residuo < tol  2.95e-14 
%   13  3.32194551e-02   2.69e+00   0.500      2  residuo < tol  1.63e-14 
%   14  1.18396588e-02   1.59e+00   1.000      1  residuo < tol  3.71e-14 
%   15  2.97211046e-03   1.17e+00   1.000      1  residuo < tol  1.27e-13 
%   16  3.60017089e-04   3.47e-01   1.000      1  residuo < tol  3.00e-14 
%   17  1.04956387e-05   8.92e-02   1.000      1  residuo < tol  8.42e-14 
%   18  1.30380006e-08   2.34e-03   1.000      1  residuo < tol  1.60e-14 
%   19  2.21713748e-14   4.21e-06   1.000      1  residuo < tol  5.00e-16 
%   20  1.33212659e-14   1.03e-07   1.000      1  residuo < tol  1.03e-07 
