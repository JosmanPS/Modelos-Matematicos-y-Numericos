function r = LBFGStlr( q, rho, s, y, H0, k, m)
 
m = min(k,m);
alpha = zeros(1,m);
index = backIndex(k,m);
 
for i = 1:m
    
    j = index(i);
    alpha(j) = rho(j) * s(:, j)' * q;
    q = q - alpha(j) * y(:, j);
    
end

r = H0 * q;
index = forwIndex(k,m);
 
for i = 1:m
    
    j = index(i);
    beta = rho(j) * y(:, j)' * r;
    r = r + s(:, j) * (alpha(j) - beta);
    
end