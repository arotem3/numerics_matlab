function [c,nf] = line_min(f,tol)
% LINE_MIN: optimization subproblem, finds local minimum on line
% f : objective function along line
% tol : relative tolerance i.e. stopping condition
% ---
% c : step size
% nf : number of function evaluations
[a,fa,b,fb,nf] = ab_init(f);
if abs(fb) < tol
    c = b;
    return
end

c = (a+b)/2; fc = f(c); nf = nf + 1;
if abs(fc) < tol
    return
end

k = 1;
while abs(fc) > tol && abs(b-a) > tol
    if k >= 100, break, end
    if (abs(fa-fc) > 1e-8) && (abs(fb-fc) > 1e-8) % inverse quaratic
        s = (a*fb*fc)/((fa-fb)*(fa-fc)) + (b*fa*fc)/((fb-fa)*(fb-fc)) + (c*fa*fb)/((fc-fa)*(fc-fb));
    else % secant
        if fc > 0
            s = c - fc*(c-a)/(fc-fa);
        else
            s = b - fb*(b-c)/(fb-fc);
        end
    end
    fs = f(s); nf = nf + 1;
    if fs < 0
        a = c;
        fa = fc;
        c = s;
        fc = fs;
    else
        b = c;
        fb = fc;
        c = s;
        fc = fs;
    end
    k = k + 1;
end
end

function [a,fa,b,fb,nf] = ab_init(f)
a = 0;
b = 1;
q = 0.75; p = 1.25;
fb = f(b);
nf = 1;
k = 1;
while fb < 0
    if k > 100, break, end
    fb1 = f(q*b); nf = nf + 1;
    if fb1 < 0 % q*b is no good
        fb2 = f(p*b); nf = nf + 1;
        if fb2 > 0 % p*b is good
            b = p*b;
            fb = fb2;
            break
        elseif fb1 > fb2 % q*b is better than p*b
            b = q*b;
            fb = fb1;
        else % p*b is better than q*b
            b = p*b;
            fb = fb2;
            [p,q] = swap(p,q); % next loop, look in direction of p first
        end
    else
        b = q*b;
        fb = fb1;
    end
    k = k + 1;
end
fa = f(a); nf = nf + 1;
end


function [y,x] = swap(x,y)
end
