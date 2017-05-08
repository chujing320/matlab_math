function [data_f, data_g, x0, feva] = plusPRP(ObjFun, x0, tol, maxiter)

    if nargin==2
        tol=1e-8;
        maxiter = 200;
    elseif nargin==3
        maxiter = 200;
    elseif nargin<2 || nargin>4
        err('error input');
    end
    
    feva = 0;beta = 0; k = 1; x1 = x0;
    [n,t] = size(x0);
    [f0, g0] = feval(ObjFun, x0, 2);
    feva = feva + 2;
    while norm(g0)>=tol
        data_f(:,k) = f0;
        daga_g(:,k) = g0;
        dk = -g0;
        [alaph, feva] = mybostwolf(ObjFun, x0, dk ,feva)
        x1 = x0 + alaph*dk;
        [f1, g1] = feval(ObjFun, x1, 2);
        feva = feva + 2;
        beta = max((g1'*(g1-g0))/(g0'*g0),0);
        dk = -g1 + beta*dk;
        k = k+1;
        if k > maxiter
            info('k>maxiter');
            break
        end
    end
end
        
        
        
    