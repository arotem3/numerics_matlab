function [x,k,X] = nlcgd(df,x,tol)
% NLCGD: nonlinear conjugate gradient method for optimization
% df : gradient function R^n -> R^n
% x : initial guess of root in R^n
% tol : relative tolerance i.e. stopping condition
% ---
% x : solution
% k : number of iterations
% X : all iterations

n = length(x);
p = -df(x);
s = p;
X = x; % -- remove

k = 0;
while norm(p,Inf) > tol
    if k >= 1000, break, end
    r = 1/norm(p,Inf);
    alpha = line_min(@(a)r*p'*df(x+a*r*p),1e-8);
    sts = s'*s;
    ds = alpha*r*p - s;
    
    s = alpha*r*p;
    if any(isnan(s)), break, end
    x = x + s;
    X = [X,x]; % -- remove
    k = k + 1;
    if mod(k,n) == 0
        p = -df(x);
    else
        beta = max([s'*ds / sts, 0]);
        p = -df(x) + beta*p;
    end
end
sp_plot(X)
end