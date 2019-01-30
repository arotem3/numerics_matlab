function r = LBFGS_update(g,S,Y,hdiag)
% LBFGS_UPDATE: determines search direction for LBFGS
% g : negative gradient at our current guess
% S : s_i = x_i - x_(i-1)
% Y : y_i = f_i - f_(i-1)
% hdiag : diagonal of hessian or estimate
% ---
% r : search direction
[~,k] = size(S);

ro = zeros(k,1);
for i=1:k
    ro(i) = 1/dot(Y(:,i),S(:,i));
end

q = g;
alpha = zeros(k,1);
beta = zeros(k,1);

for i=k:-1:1
    alpha(i) = ro(i)*dot(S(:,i),q);
    q = q - alpha(i)*Y(:,i);
end

r = hdiag.*q;

for i=1:k
    beta(i) = ro(i)*dot(Y(:,i),r);
    r = r + S(:,i)*(alpha(i) - beta(i));
end
end