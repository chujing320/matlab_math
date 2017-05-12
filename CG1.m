function [data_f, data_g, x1, k, feva] = CG1(ObjFun, x0, tol, maxiter)

    if nargin==2
        tol=1e-8;
        maxiter = 2000;
    elseif nargin==3
        maxiter = 2000;
    elseif nargin<2 || nargin>4
        err('error input');
    end
    
    rho = 0.08; sigma = 0.09; k = 1;
    [n,t] = size(x0);
    x1 = x0;feva = 0;
    g0 = ones(n,t); % 初始化g0为1矩阵
    dk = zeros(n,t); % 初始化dk为0矩阵
    [f1,g1] = feval(ObjFun, x1, 2);
    while norm(g1)>=tol
        feva = feva+2
        data_f(:,k) = f1;
        data_g(:,k) = g1;
        a = g1'*g1/(g0'*g0); %FR方法系数
        b = g1'*(g1-g0)/(g1'*g1); %PRP方法系数
        if b>a
            beta = a;
        elseif b<-a
            beta = -a;
        else
            beta = b;
        end
        if mod(k-1,50)==0 %从第0次开始算
            dk = -g1;
        else
            dk = -g1+beta*dk; 
        end
        
        f0 = f1; g0 = g1;
        [alpha,feva] = mybostwolf(ObjFun,x1,dk,feva,f0,g0,sigma,rho);
        x1 = x1+alpha*dk;
        k = k+1;
        if k > maxiter
            info('k>maxiter');
            break
        end
        [f1,g1] = feval(ObjFun, x1, 2);
    end
end