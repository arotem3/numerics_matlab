function [x,k,X,fc] = momentum_gd(f,x,tol)
% MOMENTUM_GD: momentum gradient descent
% f : gradient function R^n -> R^n
% x : initial guess of root in R^n
% tol : relative tolerance i.e. stopping condition
% ---
% x : solution
% k : number of iterations
% X : all iterations
% fc : number of function evaluations
beta = 0.8;

p = f(x);
X = x; % -- remove
r = norm(p,Inf);
[alpha,nf] = line_min( @(a) gg(a,x,-p/r,f), 1e-8);
x = x - alpha*p/r;
X = [X,x]; % -- remove
fc = 1 + nf; % num function evals

k = 1;
while norm(alpha*r,Inf) > tol
    if any(isnan(p)), break, end
    p = beta*p + f(x); fc = fc + 1;
    r = norm(p,Inf);
    [alpha,nf] = line_min( @(a) gg(a,x,-p/r,f), 1e-8);
    fc = fc + nf;
    x = x - alpha*p/r;
    X = [X,x]; % -- remove
%     if k < 30, clf(), sp_plot(X), pause(), end % -- remove
    k = k + 1;
    if k >= 1000, break, end
end
sp_plot(X) % -- remove
end

function y = gg(a,x,p,f)
z = f(x + a*p);
y = p' * z;
end