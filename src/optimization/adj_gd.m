function [x,k,X] = adj_gd(df,x,tol)
% ADJ_GD: an adjusted gradient descent hueristic.
% df : gradient function
% x : inital guess
% tol : relative tolerance i.e. stopping condition
% ---
% x : solution
% k : number of iterations
% X : all iterates
% I thought this was the same as conjugate gradient descent
% but easier to implement, turns out it isn't... that's kinda cool
% so we can call it "adjusted" gradient descent, and it works well...
X = x;
alpha = 0.1; r = 1;
k = 1;
while abs(alpha*r) > tol
    if k >= 1000, break, end % too many iterations
    prev_x = x;
    for j=1:2
        p = -df(x); r = norm(p,Inf);
        alpha = line_min(@(a) (p/r)' * df(x + (a/r)*p),1e-8);
        x = x + (alpha/r)*p;
        X = [X,x];
        k = k + 1;
    end
    p = prev_x - x; r = norm(p,Inf);
    alpha = line_min(@(a) (p/r)' * df(x + (a/r)*p),1e-8);
    x = x + (alpha/r)*p;
    X = [X,x]; % -- remove
    k = k + 1;
end
sp_plot(X)% -- remove
end