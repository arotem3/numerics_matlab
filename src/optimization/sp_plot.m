function sp_plot(X)
% SP_PLOT: spline plot for displaying iterations of root finding algorithms
% X : all iterations
[~,m] = size(X);
t = 1:m;
tt = 1:0.01:m;
xx = spline(t,X,tt);
plot(tt,xx), hold on
plot(X','o')
end