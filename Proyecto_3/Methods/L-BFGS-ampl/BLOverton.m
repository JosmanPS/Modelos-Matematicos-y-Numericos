function [ x,normg,iter, t ] = BLOverton(fun,x0,maxiter )
%Alejandra Lelo de Larrea Ibarra
%124433
%�Seguimos usando diferencias finitas? o los valores que arroja linesch_sw???

tol=10^-6;
n=length(x0);

tic;

iter=0;
x=x0;
H=eye(n);
f=feval(fun,x);
g=diffgrad(x,fun);
normg=norm(g);
curv=1;
while normg>tol && maxiter>iter && curv>0
    p=-H*g;
    [alpha, xx, ff, gg, fail, nsteps]=linesch_sw(x, f, g, p, fun, 10^-4, 0.9, 0);
    s=xx-x;
    y=gg-g;
    curv=s'*y;
    H=H+((s'*y+y'*H*y)*(s*s'))/((s'*y)^2)-((H*y*s'+s*y'*H)/(s'*y));
    normg=norm(gg,2);
    x=xx;
    g=gg;
    f=ff;
    iter=iter+1;
end

t = toc;

end

