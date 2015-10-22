function [ grad ] = diffgrad( x, fun );
%
% ... This routine computes a finite difference approximation for the
% gradient and for the Hessian matrix. The method makes use of centered 
% differences of f.
%
% j-l morales february 2015
% ITAM 
%--------------------------------------------------------------------------
%
n = length(x);
e = zeros(n,1);
grad = zeros(n,1);
hess = zeros(n,n);
eps3 = eps^(1/3);

for i=1:n
    e(i) = 1;
    if x(i) ~= 0
        h = abs(x(i))*eps3;
    else
        h = eps3;
    end
    grad(i) = (feval( fun, x + h*e ) - feval( fun, x - h*e ))/(2*h);
    e(i) = 0;
end