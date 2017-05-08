function [alpha, feva] = mybostwolfe(ObjFun, xk, dk, feva ,rho, sigma)

if nargin==4
    rho = 0.25; sigma = 0.75;
end
alpha = 1; a = 0; b = Inf; 
%while (1)
for i=0:2000
    [f0, g0] = feval(ObjFun, xk, 2);
    [f1, g1] = feval(ObjFun, xk+alpha*dk, 2);
    feva = feva +2;
    if ~(f1 <= f0 + rho*alpha*g0'*dk)
        b = alpha;
        alpha = (alpha+a)/2;
        continue;
    end
    if ~(abs(g1'*dk) <= -sigma*g0'*dk)
        a = alpha;
        alpha = min([2*alpha, (b+alpha)/2]);
        continue;
    end
    break;
end
i