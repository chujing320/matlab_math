function [alpha] = wolfe(ObjFun,xk, dk)
rho = 0.25; sigma = 0.75;
alpha = 1; a = 0; b = Inf; 
while (1)
    [f0,g0] = feval(ObjFun, xk, 2);
    [f1,g1] = feval(ObjFun, xk+alpha*dk, 2);
    if ~(f1<=f0+rho*alpha*g1'*dk)
        b = alpha;
        alpha = (alpha+a)/2;
        continue;
    end
    if ~(g1'*dk >= sigma*g0'*dk)
        a = alpha;
        alpha = min([2*alpha, (b+alpha)/2]);
        continue;
    end
    break;
end
