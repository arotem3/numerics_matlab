function J = approx_jacobian(f,x)
% APPROX_JACOBIAN: finite diff approx of jacobian of f at x.
% f : R^m -> R^n
% x : point in R^m to approx deriv at
% ---
% J : approximate jacobian [m,n]
n = length( f(x) );
m = length( x );
J = zeros(n,m);
for i=1:n
    fff = @(u) ff(u,f,i);
    J(i,:) = approx_grad(fff,x);
end
% J = sparse(J);
end

function z = ff(u,f,i)
z = f(u);
z = z(i);
end
