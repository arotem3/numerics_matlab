function g = approx_grad(f,x)
% APPROX_GRAD: finite diff approx of gradient of f at x.
% f : R^n -> R
% x : point in R^n to approx deriv at
% ---
% g : approximate gradient
n = length(x);
g = zeros(1,n);
for i=1:n
    fff = @(u) ff(u,f,x,i,n);
    g(i) = approx_deriv(fff,x(i));
end
end

function z = ff(u,f,x,i,n)
z = f([x(1:i-1); u; x(i+1:n)]);
end