function [ x, y, s ] = qp_initial_point( A, b, c, Q );

% ------------------------------------------------------
%
% Calculamos un punto inicial para el problema de 
% programación cuadrática:
%
%      min  1/2 x' * Q * x  +  c' * x
%      s.a. Ax - b = 0
%           x >= 0
%
% ------------------------------------------------------
    
    [m, n] = size(A);
    
    % Hallamos el vector de mínima norma para el problema
    % primal y dual
    
    %
    % x = A' * inv(A*A') * b;
    %
    [Q, R] = qr(A);
    x = R' \ b;
    x = Q * x;

    %
    % y = inv(A*A') * A * c;
    %
    y = Q' * c;
    y = R \ y;

    % Calculamos la holgura del problema dual
    s = Q * x + c - A' * y;
    e = ones(n, 1);

    % Hacemos positivos los vectores x y s
    delta_x = max(-3/2 * min(x), 0);
    delta_s = max(-3/2 * min(s), 0);
    x = x + delta_x * e;
    s = s + delta_s * e;
    
    % TODO (josman): Preguntar esto
    ro_x = 0.5*(x'*s)/(e'*s);
    ro_s = 0.5*(x'*s)/(e'*x);
    x  = x+ro_x*e;
    s  = s+ro_s*e;

end
