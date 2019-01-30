function [x,k] = cgd(A,b,x,tol)
% CGD: conjugate gradient method for symmetric pos def linear systems
% A*x = b <==> A'*A*x = A'*b
% A : R^m -> R^n, ideally m==n, otherwise attempts to solve lsqr soln
% x : initial guess of root in R^m
% tol : relative tolerance i.e. stopping condition
% ---
% x : solution
% k : number of iterations
if ~issymmetric(A)
    A = A'*A;
    b = A'*b;
end
r = b - A*x;
p = r;

k = 0;
while norm(r,Inf) > tol
    if k > size(A,1), break, end
    rtr = r'*r;
    alpha = rtr / (p'*A*p);
    x = x + alpha*p;
    r = r - alpha*A*p;
    beta = r'*r/rtr;
    p = r + beta*p;
    k = k + 1;
end
end