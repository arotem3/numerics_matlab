function [x,k,XX] = aafpi(g,x,m,tol)
% AAFPI: Anderson acceleration fixed point iteration.
% root finding method for problems of the form x = g(x)
% g : vector function, takes n-vec and outputs n-vec
% x : inital guess
% m : number of steps to remember
% tol : relative tolerance i.e. stopping condition
% ---
% x : solution
% k : number of iterations
% XX : all iterates
XX = x; % remove
X = []; G = [];
k = 1;
f = Inf;
while norm(f-x,Inf) > tol
    if k > 100
        break
    end
    f = g(x);
    if k <= m
        X = [X,x];
        G = [G,f];
    else
        X = [X(:,end-m+2:end), x];
        G = [G(:,end-m+2:end), f];
    end
    mk = min([k,m]);
    FF = [G - X; ones(1,mk)];
    n = size(FF,1);
    b = zeros(n,1); b(end) = 1;
    alpha = FF\b;
    
    s = G*alpha;
    x = s;
    XX = [XX,x]; % -- remove
    k = k + 1;
end
sp_plot(XX); % -- remove
end