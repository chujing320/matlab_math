function [data_f, data_g, x0, k, feva] = plusPRP(ObjFun, x0, tol, maxiter)

    if nargin==2
        tol=1e-8;
        maxiter = 2000;
    elseif nargin==3
        maxiter = 2000;
    elseif nargin<2 || nargin>4
        err('error input');
    end
    
    feva = 0;beta = 0; k = 1; x1 = x0;
    [n,t] = size(x0);
    [f0, g0] = feval(ObjFun, x0, 2);
    dk = -g0;
    data_f = 0;
    data_g = 0;
    while norm(g0)>=tol
        feva = feva + 2;
        data_f(:,k) = f0;
        daga_g(:,k) = g0;
        [alaph, feva] = mybostwolf(ObjFun, x0, dk ,feva, f0, g0, 0.5, 0.4);
        x1 = x0 + alaph*dk;
        [f1, g1] = feval(ObjFun, x1, 2);
        feva = feva + 2;
        beta = max((g1'*(g1-g0))/(g0'*g0),0);
        dk = -g1 + beta*dk;
        x0 = x1;f0 = f1;g0 = g1;
        k = k+1;
        if k > maxiter
            info('k>maxiter');
            break
        end
    end
    if feva == 0
        data_f = 0;
        data_g = 0;
        info('iter = 0');
    end
end
        
        
        
    