function df = approx_deriv(f,x)
% APPROX_DERIV: finite diff approx of deriv of f at x.
% f : one var function
% x : point to approx deriv at
% ---
% df : approximate derivative
h = 1e-2; err = 1e-6; tol = 1e-4;
df = (f(x - 2*h) - 8*f(x-h) + 8*f(x+h) - f(x+2*h))/(12*h);
if abs(df) < tol % measure of sparsity
    df = 0;
    return
end
h = h*0.75;
df1 = (f(x - 2*h) - 8*f(x-h) + 8*f(x+h) - f(x+2*h))/(12*h);
while abs(df1-df) > err
    df = df1;
    h = h*0.75;
    df1 = (f(x - 2*h) - 8*f(x-h) + 8*f(x+h) - f(x+2*h))/(12*h);
end
end