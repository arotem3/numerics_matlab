function [x,k,X] = leven_marq(f,x,tol)
% LEVEN_MARQ: trust region method for nonlinear least squares
% f : R^m -> R^n
% x : initial guess of root in R^m
% tol : relative tolerance i.e. stopping condition
% ---
% x : solution
% k : number of iterations
% X : all iterations
n = length(x);
nu = 2;
tau = 1e-3;
F = f(x);%delta = 1;
X = x; % -- remove
J = approx_jacobian(f,x);
lam = tau*max(diag(J'*J));
k = 0;
while norm(F,Inf) > tol
    if k >= 100, break, end
    RHS = -J'*F;
    rho = -1;
    while rho < 0
        A = J'*J + lam*eye(n);
        delta = A\RHS;
        F1 = f(x + delta);
        rho = (norm(F) - norm(F1))/(delta'*(lam*delta+RHS));
        if rho > 0
            x = x+delta;
            X = [X,x]; % -- remove
            J = approx_jacobian(f,x);
            F = F1;
            lam = lam * max(0.33, 1 - (2*rho-1)^3);
            nu = 2;
        else
            lam = lam*nu;
            nu = 2*nu;
        end
    end
    k = k + 1;
end
sp_plot(X)
end