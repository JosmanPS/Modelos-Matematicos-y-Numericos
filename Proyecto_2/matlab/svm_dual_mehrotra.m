function [ x, y, s ] = svm_dual_mehrotra(X, y, c, tol)

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

    varlist = {'X', 'y', 'Y', 'c'};
    clear(varlist{:});
    
    %
    % Calculamos el punto inicial para el mpi
    %
    e_n = ones(n, 1);
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
    S = diag(s);           S = sparse(S);
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
    F = - [F1; F3; -F2];

    %
    % Comenzamos el proceso iterativo
    %
    while d_gap > tol & iter < 20
        
        iter = iter + 1

        %
        % Planteamos el sistema predictor
        %
        KKT = [  zeros(n,n),       A'    ,     -b      ;
                     A     ,    -eye(m)  ,  zeros(m, 1);
                    -b'    ,  zeros(1, m),     0.0     ];
        KKT = sparse(KKT);
        cond(KKT)

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
        mu_aff = (x + alpha_x*dx)' * (s + alpha_s*ds) / n;
        sigma = (mu_aff/mu)^3;

        % % % % % % % % % % % % % % % %
        
        %
        % Planteamos el sistema corrector
        %
        z = mu * sigma * X_inv * e_n;
        w = mu * sigma * S_inv * e_n;
        Z = diag(z);    Z = sparse(Z);
        W = diag(w);    W = sparse(W);

        KKT(1:n, 1:n) = X_inv*Z + S_inv*W;
        F5 = X*z - mu * sigma * e_n;
        F6 = S*w - mu * sigma * e_n;
        F(1:n) = - (F1 - z + w + X_inv*F5 + S_inv*F6 + S_inv*W*F4);

        %
        % Resolvemos el sistema corrector
        %
        d = backsolve(KKT, F);
        dx = d(1:n);
        dy = d(n+1:n+m);
        dlm = d(n+m+1);
        ds = -F4 - dx;

        alpha_x = step(x, dx, 1);
        alpha_s = step(s, ds, 1);
        alpha = min(alpha_x, alpha_s)
        
        %
        % Calculamos el paso
        %
        x = x + alpha*dx;
        y = y + alpha*dy;
        s = s + alpha*ds;
        lm = lm + alpha*dlm;

        %
        % Actualizamos los nuevos valores
        %
        %
        X = sparse(diag(x));
        S = sparse(diag(s));
        X_inv = sparse(diag(1./x));
        S_inv = sparse(diag(1./s));

        F1 = -e_n - lm*b + A'*y;
        F2 = b' * x;
        F3 = A*x - y;
        F4 = x - gamma*e_n + s;
        F = - [F1; F3; -F2];

        mu = (x' * s) / n;
        d_gap = mu

    end
    
end
