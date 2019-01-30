% A = [9.9 5.9 7.7;...
%      5.9 11.7 5.9;...
%      7.7 5.9 9.9];
% A = [2 -1 0;...
%      -1 2 -1;...
%      0 -1 2];
A = [2.58, 2.18, 2.32;...
     2.18, 2.28, 1.56;...
     2.32, 1.56, 3.44];
f = @(x) 0.5 * x' * A * x;
df = @(x) A*x;

m = 20;
x = [1;2;-1];
X = x;
for i=1:m
    p = -df(x);
    g = @(a) gg(a,x,p,df);
    alpha = fzero(g,1e-2);
    x = x + alpha*p;
    X = [X,x];
end
plot3(X(1,:),X(2,:),X(3,:),'-or'), hold on
xlim([-2,2]), ylim([-2,2]), zlim([-2,2])

[V,D] = eig(A);
t = linspace(-10,10);
v1 = V(:,1)*t;
v2 = V(:,2)*t;
v3 = V(:,3)*t;
plot3(v1(1,:), v1(2,:),v1(3,:),'r')
plot3(v2(1,:), v2(2,:),v2(3,:),'k')
plot3(v3(1,:), v3(2,:),v3(3,:),'m')


[x,k,X] = adj_gd(df,[1;2;-1],1e-4);
x,k
plot3(X(1,:),X(2,:),X(3,:),'-ob')

figure
plot((V'*X)','-o')
legend(num2str(D(1,1)),num2str(D(2,2)),num2str(D(3,3)))

function y = gg(a,x,p,df)
y = p' * df(x + a*p);
end