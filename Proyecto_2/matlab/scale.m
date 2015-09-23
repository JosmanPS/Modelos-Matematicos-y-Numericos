function [ X ] = scale(X);
[ m, n ] = size(X); e = ones(m,1);

for i=1:n
   xm = mean( X(:,i) );
   xs = std( X(:,i) );
   X(:,i) = (X(:,i) - xm*e)/xs;
end