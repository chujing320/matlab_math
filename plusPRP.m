function [data_f, data_g, x0, k, feva] = plusPRP(ObjFun, x0, tol, maxiter)

    if nargin==2
        tol=1e-20;
        maxiter = 2000;
    elseif nargin==3
        maxiter = 2000;
    elseif nargin<2 || nargin>4
        err('error input');
    end
    
    feva = 0;beta = 0; k = 1; x1 = x0;f0=0;
    rho = 0.15; sigma = 0.35;
    [f1, g1] = feval(ObjFun, x0, 2);
  %  while norm(g1)>=tol
    while abs(f1-f0)>=tol
        feva = feva + 2;
        data_f(:,k) = f1;
        data_g(:,k) = g1;
        f0 = f1;g0=g1;x0=x1;
        if mod(k-1,50)==0 %从第0次开始算
            dk = -g1;
        else
            dk = -g1+beta*dk; 
        end
        [alaph, feva] = mybostwolf(ObjFun, x0, dk ,feva, f0, g0, sigma, rho);
        x1 = x0 + alaph*dk;
        [f1, g1] = feval(ObjFun, x1, 2);
        feva = feva + 2;
        beta = max((g1'*(g1-g0))/(g0'*g0),0);
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
        
        
        
    