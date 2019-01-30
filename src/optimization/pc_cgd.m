function [x,k] = pc_cgd(A,Mi,b,x,tol)
% PC_CGD: preconditioned conjugate gradient method for
% symmetric pos def linear systems
% A*x = b <==> A'*A*x = A'*b
% A : R^m -> R^n, ideally m==n, otherwise attempts to solve lsqr soln
% Mi : inverse of preconditioning matrix
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
z = Mi*r;
p = z;

k = 0;
while norm(r,Inf) > tol
    if k > size(A,1), break, end
    alpha = r'*z / (p'*A*p);
    x = x + alpha*p;
    ztr = z'*r;
    r = r - alpha*A*p;
    z = Mi * r;
    beta = z'*r / ztr;
    p = z + beta*p;
    k = k + 1;
end