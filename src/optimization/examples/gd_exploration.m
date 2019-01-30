% A = [3.4 2.2;...
%      2.2 2.6];
A = [1,0;0,3];
f = @(x) 0.5 * x' * A * x;
df = @(x) A*x;

% b = 100;
% f = @(x) (1-x(1))^2 + b*(x(2) - x(1)^2)^2;
% df = @(x) [-2*(1-x(1)) - 4*b*x(1)*(x(2)-x(1)^2);...
%            2*b*(x(2)-x(1)^2)];
       
% f = @(x) sum(polyval([0.75,3,0,-12,0],x));
% df = @(x) polyval([3,9,0,-12],x);

% f = @(x) sum(x.^4 - 16*x.^2 + 5*x)/2;
% df = @(x) 2*x.^3 - 16*x + 2.5;
% A = approx_jacobian(df,[1;1]);
        
a = -5; c = 5; b = -5; d = 5; rx = 0; ry = rx;
n = 100;
xx = linspace(a,c,n); yy = linspace(b,d,n);
[xx,yy] = meshgrid(xx,yy);
zz = zeros(size(xx));
for i=1:n
    for j=1:n
        zz(i,j) = f([xx(i,j);yy(i,j)]);
    end
end

contourf(xx,yy,zz,[0:0.002,1,1:2:f([5;5])]), hold on
xlim([a c]), ylim([b d])

[V,~] = eig(A);
t = linspace(-10,10);
v1 = V(:,1)*t + rx;
v2 = V(:,2)*t + ry;
plot(v1(1,:), v1(2,:),'r')
plot(v2(1,:), v2(2,:),'r')

m = 200;
[x,y] = ginput(1); x0 = [x;y]; x = x0;

B = approx_jacobian(df,x0);
[V,~] = eig(B);
t = linspace(-10,10);
v1 = V(:,1)*t + x0;
v2 = V(:,2)*t + x0;
plot(v1(1,:), v1(2,:),'g')
plot(v2(1,:), v2(2,:),'g')

X = x;
for i=1:m
    p = -df(x);
    g = @(a) gg(a,x,p,df);
    alpha = fzero(g,1e-2);
    x = x + alpha*p;
    X = [X,x];
end
plot(X(1,:),X(2,:),'-or')

x = x0;
[~,~,X] = adj_gd(df,x,1e-4);
plot(X(1,:),X(2,:),'-*y')

% x = x0;
% X = x;
% for i=1:m/2
%     p = -df(x);
%     p = [p(1);0];
%     g = @(a) gg(a,x,p,df);
%     alpha = fzero(g,1e-2);
%     x = x + alpha*p;
%     X = [X,x];
%     
%     p = -df(x);
%     p = [0;p(2)];
%     g = @(a) gg(a,x,p,df);
%     alpha = fzero(g,1e-2);
%     x = x + alpha*p;
%     X = [X,x];
% end
% plot(X(1,:),X(2,:),'-sw')
% 
% x = x0;
% X = x;
% p = -df(x);
% p = [p(1);0];
% g = @(a) gg(a,x,p,df);
% alpha = fzero(g,1e-2);
% x = x + alpha*p;
% X = [X,x];
% 
% p = -df(x);
% p = [0;p(2)];
% g = @(a) gg(a,x,p,df);
% alpha = fzero(g,1e-2);
% x = x + alpha*p;
% X = [X,x];
% for i=3:m/2
%     for j=1:2
%         p = -df(x);
%         g = @(a) gg(a,x,p,df);
%         alpha = fzero(g,1e-2);
%         x = x + alpha*p;
%         X = [X,x];
%     end
%     p = X(:,i) - x;
%     g = @(a) gg(a,x,p,df);
%     alpha = fzero(g,1e-2);
%     x = x + alpha*p;
%     X = [X,x];
% end
% plot(X(1,:),X(2,:),'-dc')

function y = gg(a,x,p,df)
y = p' * df(x + a*p);
end