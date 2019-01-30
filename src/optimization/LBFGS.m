function [x,k,X] = LBFGS(obj_func,f,x0,m,tol)
% LBFGS: limited memory BFGS Quasi newton method for optimization
% obj_func : objective function
% f : gradient function R^n -> R^n
% x : initial guess of root in R^n
% m : number of steps to remember
% tol : relative tolerance i.e. stopping condition
% ---
% x : solution
% k : number of iterations
% X : all iterations
k = 1;
x = x0;
X = x; % -- remove

n = length(x0);
hdiag = ones(n,1);
p = -f(x);
alpha = wolfe_cond(x,p,obj_func,f);
% alpha = fzero( @(a) p' * f(x + a*p), 0);
s = alpha*p;
x = x + s;
X = [X,x]; % -- remove
y = f(x) - f(x0);

S = s;
Y = y;

while norm(p,Inf) > tol
    p = LBFGS_update(-f(x),S,Y,hdiag);
    alpha = wolfe_cond(x,p,obj_func,f);
%     alpha = fzero( @(a) p' * f(x + a*p), 0);
    s = alpha*p;
    x = x + s;
    X = [X,x]; % -- remove
%     clf(), sp_plot(X), pause() % -- remove
    y = f(x) - f(x-s);
    if k < m+1 % instead we should push onto cyclic queue
        S = [S,s];
        Y = [Y,y];
    else
        S = [S(:,2:m-1),s];
        Y = [Y(:,2:m-1),y];
    end
    k = k + 1;
end
sp_plot(X)
end