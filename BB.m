function [data_f, data_g, x0, feva] = BB(ObjFun, x0, tol, maxiter)

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
    feva = feva + 1;
    alaph = xxxxxxxxxxxx;
    x1 = x0 - alaph*g0;
    data_f(:,k) = f0;
    daga_g(:,k) = g0;
    k = 2;
    while norm(g0)>=tol
        [f1, g1] = feval(ObjFun, x1, 2);
        feva = feva + 1;
        data_f(:,k) = f1;
        daga_g(:,k) = g1;
        sk = x1-x0; yk = g1-g0;
        alaph = (sk'*sk)/(sk'*yk);
        x0 = x1; g0 = g1;
        x1 = x1 - alaph*g1;    
        k = k+1;
        if k > maxiter
            info('k>maxiter');
            break
        end
    end
end
        
        
        
    