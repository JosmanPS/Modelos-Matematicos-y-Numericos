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
    lm = 1;
    z = x;
    w = x;
    mu = (x'*z + s'*w) / (2*n);
    d_gap = mu;

    %
    % Ajustamos los otros parámetros necesarios para el modelo
    %
    X = diag(x);           X = sparse(X);
    S = diag(s);           S = sparse(S);
    X_inv = diag(1./x);    X_inv = sparse(X_inv);
    S_inv = diag(1./s);    S_inv = sparse(S_inv);
    Z = diag(z);           Z = sparse(Z);
    W = diag(w);           W = sparse(W);
    iter = 0;

   
    % 
    % Escribimos los residuos sin parámetro de centralidad
    %
    F1 = -e_n - lm*b + A'*y -z + w;
    F2 = b' * x;
    F3 = A*x - y;
    F4 = x - gamma*e_n + s;
    F5 = X*z;
    F6 = S*w;
    F1_cor = F1 + X_inv*F5 - S_inv*F6 + S_inv*W*F4;
    F = - [F1+X_inv*F5-S_inv*F6+S_inv*W*F4; F3; -F2];
    obj = 0.5 * y'*y - e_n'*x;

    %
    % Impresiones iniciales
    %
    fprintf(['iter   ||f1||      ||f2||      ||f3||      ||f4||      ||f5||      ||f6||   ' ...
             '     OBJ        mu_k        alpha       sigma\n']);
    fprintf(['---------------------------------------------------------------' ...
             '---------------------------------------------------------------\n']);
    
    %
    % Comenzamos el proceso iterativo
    %
    while d_gap > tol && iter < 50 && norm(F) > tol
        
        iter = iter + 1;

        %
        % Planteamos el sistema predictor
        %
        KKT = [  X_inv*Z + S_inv*W,       A'    ,     -b      ;
                         A        ,    -eye(m)  ,  zeros(m, 1);
                        -b'       ,  zeros(1, m),     0.0     ];
        KKT = sparse(KKT);

        %
        % Resolvemos el sistema predictor
        %
        d = backsolve(KKT, F);
        dx = d(1:n);
        dy = d((n+1):(n+m));
        ds = -F4 - dx;
        dlm = d(n+m+1);
        dz = -X_inv * (F5 + Z*dx);
        dw = -S_inv * (F6 + W*ds);

        alpha_x = step(x, dx, 1);
        alpha_s = step(s, ds, 1);
        alpha_z = step(z, dz, 1);
        alpha_w = step(w, dw, 1);

        %
        % Calculamos el parámetro de centralidad
        %
        mu_aff = (x + alpha_x*dx)' * (z + alpha_z*dz);
        mu_aff = mu_aff + (s + alpha_s*ds)' * (w + alpha_w*dw);
        mu_aff = mu_aff / (2*n);
        sigma = (mu_aff/mu)^3;

        % % % % % % % % % % % % % % % %
        
        %
        % Planteamos el sistema corrector
        %
        % z = mu * sigma * X_inv * e_n;
        % w = mu * sigma * S_inv * e_n;
        % Z = diag(z);    Z = sparse(Z);
        % W = diag(w);    W = sparse(W);

        % KKT(1:n, 1:n) = X_inv*Z + S_inv*W;
        F5 = F5 - mu * sigma * e_n;
        F6 = F6 - mu * sigma * e_n;
        F1_cor = F1 + X_inv*F5 - S_inv*F6 + S_inv*W*F4;
        F = - [F1_cor; F3; -F2];

        %
        % Resolvemos el sistema corrector
        %
        d = backsolve(KKT, F);
        dx = d(1:n);
        dy = d(n+1:n+m);
        dlm = d(n+m+1);
        ds = -F4 - dx;
        dz = -X_inv * (F5 + Z*dx);  %%%%
        dw = -S_inv * (F6 + W*ds);  %%%%

        %
        % Calculamos el paso
        %
        alpha_x = step(x, dx, 0.9995);
        alpha_s = step(s, ds, 0.9995);
        alpha_z = step(z, dz, 0.9995);
        alpha_w = step(w, dw, 0.9995);
        alpha = min(alpha_x, alpha_s);
        
        %
        % Actualizamos los vectores
        %
        x = x + alpha_x*dx;
        y = y + alpha*dy;
        s = s + alpha_s*ds;
        lm = lm + alpha*dlm;
        z = z + alpha_z*dz;
        w = w + alpha_w*dw;

        %
        % Actualizamos los nuevos valores
        %
        X = sparse(diag(x));
        S = sparse(diag(s));
        X_inv = sparse(diag(1./x));
        S_inv = sparse(diag(1./s));
        Z = diag(z);    Z = sparse(Z);
        W = diag(w);    W = sparse(W);

        F1 = -e_n - lm*b + A'*y -z + w;
        F2 = b' * x;
        F3 = A*x - y;
        F4 = x - gamma*e_n + s;
        F5 = X*z;
        F6 = S*w;
        F1_cor = F1 + X_inv*F5 - S_inv*F6 + S_inv*W*F4;
        F = - [F1+X_inv*F5-S_inv*F6+S_inv*W*F4; F3; -F2];

        mu = (x'*z + s'*w) / (2*n);
        d_gap = mu;
        obj = 0.5 * y'*y - e_n'*x;

        fprintf(['%3i  %1.4e  %1.4e  %1.4e  %1.4e  %1.4e  %1.4e  %1.4e  %1.4e  ' ...
                 '%1.4e  %1.4e \n'], iter, norm(F1), norm(F2), norm(F3), norm(F4), ...
                norm(F5), norm(F6), obj, mu, alpha, sigma);

    end
    
end