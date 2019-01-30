function conjugate_grad_compare(f,x0,edges)
% exploration into nonlinear optimization comparing conjugate gradient
% descent and my adjusted gradient descent. 
if nargin == 0
    f = @df;
    a = -2; b = 2;
    edges = -2:0.1:2;
    x0 = (b-a)*rand(n,1) + a;
end
if nargin == 2
    edges = 100;
end
tol = 1e-5;


subplot(2,2,1)
disp("using nlcgd()...")
tic
[x,k] = nlcgd(f,x0,tol);
toc
histogram(x,edges)
subplot(2,2,3)
plot(x,'-o')
disp("k = " + num2str(k) + "    x = " + num2str(mean(x)) + "    ||df|| = " + norm(f(x),Inf) );

subplot(2,2,2)
disp("using adj_gd()...")
tic
[x,k] = adj_gd(f,x0,tol);
toc
histogram(x,edges)
subplot(2,2,4)
plot(x,'-o')
disp("k = " + num2str(k) + "    x = " + num2str(mean(x)) + "    ||df|| = " + norm(f(x),Inf) );
end

function y = df(x)
D = length(x);
y = zeros(D, 1);
y(1:D-1) = - 400*x(1:D-1).*(x(2:D)-x(1:D-1).^2) - 2*(1-x(1:D-1));
y(2:D) = y(2:D) + 200*(x(2:D)-x(1:D-1).^2);
end