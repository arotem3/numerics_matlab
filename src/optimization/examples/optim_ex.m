clc, close all
str = "(1) rosenbrock" + newline + "(2) styblinski tang" + newline + "(3) polynomial" + newline + "please make a selection: ";
i = input(char(str));
figure(1)
set(gcf,'position',[100,100,800,600])
if i == 1
    % -- rosenbrock
    disp(newline + "f(x,y) = (1-x)^2 + 100(y - x^2)^2" + newline)
    b = 100;
    c = 1; n = 2;
    xx = linspace(-2,2,1000); yy = linspace(-1,3,1000);
        [xx,yy] = meshgrid(xx,yy);
        zz = (1-xx).^2 + b*(yy - xx.^2).^2;
        contourf(xx,yy,zz,[0,2.^([0:7,7:0.5:11])]), hold on
        plot(1,1,'or','MarkerSize',20)
        xlim([-2,2]), ylim([-1,3])
    [x,y] = ginput(1); x0 = [x;y];
    f = @(x) (1-x(1))^2 + b*(x(2) - x(1)^2)^2;
    df = @(x) [-2*(1-x(1)) - 4*b*x(1)*(x(2)-x(1)^2);...
              2*b*(x(2)-x(1)^2)];
elseif i == 2
    % -- styblinsky tang
    disp(newline + "f(x,y) = sum(x^4 - 16x^2 + 5x)/2" + newline)
    xx = linspace(-5,5);
        [xx,yy] = meshgrid(xx);
        zz = (xx.^4 - 16*xx.^2 + 5*xx)/2 + (yy.^4 - 16*yy.^2 + 5*yy)/2;
        contourf(xx,yy,zz,30), hold on
        plot(-2.9035,-2.9035,'or','MarkerSize',20)
        xlim([-5,5]), ylim([-5,5])
    n = 2;
    f = @(x) sum(x.^4 - 16*x.^2 + 5*x)/2;
    df = @(x) 2*x.^3 - 16*x + 2.5;
    [x,y] = ginput(1); x0 = [x;y];
elseif i == 3
    % -- polynomial with saddle point as x = -2
    disp(newline + "f(x,y) = sum(0.75*x^4 + 3x^3 - 12x)" + newline)
    n = 2;
    p = [0.75,3,0,-12,0];
    dp = [3,9,0,-12];
    xx = linspace(-5,4);
        [xx,yy] = meshgrid(xx);
        zz = polyval(p,xx) + polyval(p,yy);
        contourf(xx,yy,zz,30), hold on
        plot(1,1,'or','MarkerSize',25,'LineWidth',2)
        xlim([-5,4]), ylim([-5,4])
    f = @(x) sum(polyval(p,x));
    df = @(x) polyval(dp,x);
    [x,y] = ginput(1); x0 = [x;y];
else
    close all
    optim_ex
end
 
% -- solver parameters
steps_to_remem = 3;
use_diag = false;
tol = 1e-5;

% -- the meat
figure(2)
disp("nlcgd")
subplot(2,2,1)
tic, [x,k,X] = adj_gd(df,x0,tol); toc
% tic, [x,k,X] = nlcgd(df,x0,tol); toc
title("nlcgd")
figure(1), plot(X(1,:),X(2,:),'-ow'), figure(2)
disp("k = " + num2str(k) + "    x = " + num2str(mode(x)) + "    ||df|| = " + num2str(norm(df(x),Inf)))

disp(newline + "LBFGS")
figure(2),subplot(2,2,2)
tic, [x,k,X] = LBFGS(f,df,x0,steps_to_remem,tol); toc
title("LBFGS")
figure(1), plot(X(1,:),X(2,:),'-or','MarkerSize',10), title("Quasi Newton Optimization"), figure(2)
disp("k = " + num2str(k) + "    x = " + num2str(mode(x)) + "    ||df|| = " + num2str(norm(df(x),Inf)))

disp(newline + "leven_marq")
subplot(2,2,3)
tic, [x,k,X] = leven_marq(df,x0,tol); toc
title("leven marq")
figure(1), plot(X(1,:),X(2,:),'-oy'), figure(2)
disp("k = " + num2str(k) + "    x = " + num2str(mode(x)) + "    ||df|| = " + num2str(norm(df(x),Inf)))

disp(newline + "momentum_gd")
subplot(2,2,4)
tic, [x,k,X,fc] = momentum_gd(df,x0,tol); toc
title('momentum gd')
figure(1), plot(X(1,:),X(2,:),'-om'), figure(2)
disp("k = " + num2str(fc) + "    x = " + num2str(mode(x)) + "    ||df|| = " + num2str(norm(df(x),Inf)))

figure(1)
legend({"objective function","minima","nlcgd","LBFGS","leven marq","momentum gd"},'Location','bestoutside','Color','k','TextColor','w')

function d = ff(df,x)
m = length(x); n = m/2 + 1;
a = -df(x(1:n-1)) - 0.99*x(n:end);
d = [x(n:end);a];
end

function [val,stop,dir] = ev(t,x,tol)
m = length(x)/2 + 1;
val = norm( x(m:end),Inf ) > tol; val = val || (t<tol); val = val*1;
stop = 1;
dir = 0;
end