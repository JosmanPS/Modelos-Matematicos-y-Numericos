function [ x, y, s ] = ipm(X, y, c)

    %
    % Calculamos los valores del modelo computacional
    %
    %      min    1/2 y' y  -  e' x
    %      s.a.   b' x = 0
    %             A x - y = 0
    %             x - gamma*e + s = 0
    %             x >= 0
    %             s >= 0

    Y = diag(y);
    [n_samples, n_features] = size(X)
    X = X';
    A = X * Y;
    b = y;
    gamma = c;

    clear(X, y, Y, c);
    
    %
    % Calculamos el punto inicial para el mpi
    %

    [ x, y, s ] = qp_inital_point(A, b, c , Q);
    