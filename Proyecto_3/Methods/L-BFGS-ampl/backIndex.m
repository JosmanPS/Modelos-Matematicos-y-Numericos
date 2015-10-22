function index = backIndex(k, m)

if k < m
    
    index = (k-1):-1:1;
    
else

    index = zeros(m,1);
    n = mod(k-1, m);
    index(1:n) = [n:-1:1];
    index((n+1):m) = [m:-1:n+1];

end
    
end

