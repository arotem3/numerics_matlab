function alpha = wolfe_cond(x,p,obj_func,f)
% WOLFE_COND: approximately solves line min problem using wolfe condition
% x : current guess of root in R^n
% p : search direction in R^n
% obj_func : objective function
% f : gradient function R^n -> R^n
% ---
% alpha : step size
wolfe_c1 = 1e-4;
wolfe_c2 = 0.9;
b = 0.5;
phi = @(a) obj_func(x + a.*p);
dphi = @(a) f(x + a.*p);

alpha = 1;

cc = 0;
while true
    pdfx = p'*dphi(alpha);
    cond1 = phi(alpha) <= obj_func(x) + wolfe_c1*alpha*pdfx;
    cond2 = abs(p'*dphi(alpha)) <= wolfe_c2*abs(p'*dphi(0));
    if cond1 && cond2
        break;
    elseif phi(b*alpha) < phi(alpha)
        alpha = (1-b)*alpha;
    else
        alpha = (1+b)*alpha;
    end
    cc = cc + 1;
    if cc >= 1000, break, end % stoping criteria
end
end