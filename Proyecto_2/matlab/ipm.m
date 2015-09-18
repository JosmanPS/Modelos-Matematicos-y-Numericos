function [ x, y, s ] = ipm(X, y, c=1, tol=1e-6)

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
    Y = sparse(Y);    
    X = X';
    A = X * Y;
    b = y;
    gamma = c;
    [m, n] = size(A);

    clear(X, y, Y, c);
    
    %
    % Calculamos el punto inicial para el mpi
    %
    e_n = zeros(n, 1);
    x = gamma / 2 * e_n;
    s = x;
    y = A * x;
    mu = (x' * s) / n;
    d_gap = mu;
    lm = 1;

    %
    % Ajustamos los otros parámetros necesarios para el modelo
    %
    X = diag(x);           X = sparse(X);
    Y = diag(y);           Y = sparse(Y);
    X_inv = diag(1./x);    X_inv = sparse(X_inv);
    S_inv = diag(1./s);    S_inv = sparse(S_inv);
    iter = 0;

    % 
    % Escribimos los residuos sin parámetro de centralidad
    %
    F1 = -e_n - lm*b + A'*y;
    F2 = b' * x;
    F3 = A*x - y;
    F4 = x - gamma*e_n + s;

    %
    % Comenzamos el proceso iterativo
    %
    while d_gap > tol & iter < 20
        
        iter = iter + 1;

        %
        % Planteamos el sistema a resolver
        %
        KKT = [  zeros(n,n),       A'    ,     -b      ;
                     A     ,    -eye(m)  ,  zeros(m, 1);
                    -b'    ,  zeros(1, m),     0.0     ];
        KKT = sparse(KKT);
        F = - [F1; F3; -F2];

        %
        % Resolvemos el sistema predictor
        %
        d = backsolve(KKT, F);
        dx = d(1:n);
        ds = -F4 - dx;

        alpha_x = step(x, dx, 1);
        alpha_s = step(s, ds, 1);

        %
        % Calculamos el parámetro de centralidad
        %
        