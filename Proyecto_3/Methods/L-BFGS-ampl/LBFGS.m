function x = LBFGS(fun, maxiter, m, tol, mod_print, graphics)

tic
nombre_ampl = strcat(fun, '.nl');

[x, ~, ~, ~, ~, ~] = spamfunc(nombre_ampl);

% Valores iniciales
n = length(x);
k = 1;
H = 1;
[f,~] = spamfunc(x,0);
[g,~] = spamfunc(x,1);
normg = norm(g);
curv = 1;
S = zeros(n,m);
Y = zeros(n,m);
R = zeros(m,1);

% Preparamos la impresion
fprintf(1,'   k        f            ||g||     alfa     nfg     curv  \n');
fprintf(1,'-----------------------------------------------------------------------\n');
fprintf(' %3i  %14.8e   %8.2e  \n', ...
                            k, f,  normg );

% Almacenamos los primeros m datos
%while k < m + 1 && normg > tol && curv > 0
%    
p= -H * g;
[alpha, xx, ff, gg, ~, numfg] = linesch_swAMPL(x, f, g, p, fun, 10^-4, 0.9, 0);
s = xx - x;
y = gg - g;
curv = s' * y;
  %  H = H +((s'*y+y'*H*y)*(s*s'))/((s'*y)^2)-((H*y*s'+s*y'*H)/(s'*y));
normg = norm(gg,Inf);
x = xx;
g = gg;
f = ff;
S(:, k) = s;
Y(:, k) = y;
R(k) = 1 / (y' * s);
k = k + 1;
fprintf(' %3i  %14.8e   %8.2e   %5.3f    %3i       %5.3f   \n', ...
                            k, f,  normg, alpha, numfg, curv );

%end
graph = [];
H = s'*y/(y'*y);
% Usamos la implementacion de memoria limitada
while k < maxiter && normg > tol && curv > 0
    
    graph = [graph, normg];
    
    p = LBFGStlr(g, R, S, Y, H, k, m);
    p = -p;
    [alpha, xx, ff, gg, ~, numfg] = linesch_swAMPL(x, f, g, p, fun, 10^-4, 0.9, 0);
    s = xx - x;
    y = gg - g;
    curv = s' * y;
    normg = norm(gg, Inf);
    x = xx;
    g = gg;
    f = ff;
    
    ind = mod(k,m);
    if ind == 0
        ind = m;
    end
    
    S(:, ind) = s;
    Y(:, ind) = y;
    R(ind) = 1 / (y' * s);
    k = k + 1;
    if mod(k, mod_print) == 0
        fprintf(' %3i  %14.8e   %8.2e   %5.3f    %3i       %5.3f   \n', ...
                            k, f,  normg, alpha, numfg, curv );
    end
    
    
end
fprintf(' %3i  %14.8e   %8.2e   %5.3f    %3i       %5.3f   \n', ...
                            k, f,  normg, alpha, numfg, curv );
toc

if graphics
    semilogy(graph)
end

end