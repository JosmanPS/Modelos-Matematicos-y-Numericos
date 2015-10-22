function index = forwIndex(k, m)

if k < m
    
    index = 1:(k-1);
    
else

    index = zeros(m,1);
    n = mod(k, m);
    if n == 0
        n = m;
    end
    index(1:m-n+1) = [n:m];
    index((m-n+2):m) = [1:n-1];

end
    
end
