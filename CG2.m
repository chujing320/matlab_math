function [data_f, data_g, x1, k, feva] = CG2(ObjFun, x0, tol, maxiter)

    if nargin==2
        tol=1e-8;
        maxiter = 2000;
    elseif nargin==3
        maxiter = 2000;
    elseif nargin<2 || nargin>4
        err('error input');
    end
    
    rho = 0.18; sigma = 0.68; k = 1;
    [n,t] = size(x0);
    x1 = x0; feva =0;
    g0 = ones(n,t); % 初始化g0为1矩阵
    dk = zeros(n,t); % 初始化dk为0矩阵
    [f1,g1] = feval(ObjFun, x1, 2);
    while norm(g1)>=tol
        feva = feva+2;
        data_f(:,k) = f1;
        data_g(:,k) = g1;
        a = g1'*g1/(g0'*g0); %FR方法系数
        b = g1'*(g1-g0)/(g1'*g1); %PRP方法系数
    %    if mod(k-1,50)==0 %从第0次开始算 n步重新开始准则
     %       dk = -g1;
     %   else
            dk = -g1+max(0,min(a,b))*dk;
      %  end  
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
    if feva == 0
        data_f = 0;
        data_g = 0;
        info('iter = 0');
    end
end