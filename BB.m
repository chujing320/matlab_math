 function [data_f, data_g, x0, k, feva] = BB(ObjFun, x0, tol, maxiter)

    if nargin==2
        tol=1e-14;
        maxiter = 2000;
    elseif nargin==3
        maxiter = 2000;
    elseif nargin<2 || nargin>4
        err('error input');
    end
    
    feva = 0; k = 1; 
    [n,t] = size(x0);
    [f1, g1] = feval(ObjFun, x0, 2);   
    x1 = x0;
    x0 = zeros(n,t);
    g0 = zeros(n,t);
    while norm(g1)>=tol
        feva = feva + 2;
        data_f(:,k) = f1;
        data_g(:,k) = g1;
        sk = x1-x0; yk = g1-g0;
        alaph = (sk'*sk)/(sk'*yk+10^-12);
        x0 = x1; g0 = g1;
        x1 = x1 - alaph*g1; 
        k = k+1;
        if k > maxiter
            info('k>maxiter');
            break
        end
        [f1, g1] = feval(ObjFun, x1, 2);
    end
    if feva == 0
        data_f = 0;
        data_g = 0;
        info('iter = 0');
    end
end
        
        
        
    