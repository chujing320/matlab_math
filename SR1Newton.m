function [data_f,data_g, x0,feva] = SR1Newton(ObjFun,x0,tol,maxiter)

    if nargin==2
        tol=1e-8;
        maxiter = 200;
    elseif nargin==3
        maxiter = 200;
    elseif nargin<2 || nargin>4
        err('error input');
    end
    [n,t] = size(x0);

    k=1;
    x1=x0;
    H0 = ones(n,n);
    H1 = H0
    [f0 g0]=feval(ObjFun, x0, 2);
    feva = 2;
    while norm(g0)>=tol
        %data 用来存放中间数据
        data_f(:,k) = f0;
        data_g (:,k) = g0;
        %开始迭代
        d = -H0*g0; 
        [alaph,info1] = bolinesearch(ObjFun, x0, d, Rule);
        x1 = x0+alaph*d;%线搜索准则
        [f0 g0]=feval(ObjFun, x1, 2);%传入返回值个数2
        feva = feva+3;
        %修正公式
        s = x1-x0;
        y = g1-g0;
        H1 = (H0+(s-H0*y)*(s-H0*y)')/((s-H0*y)'*y);
        x0 = x1;
        g0 = g1;
        H0 = H1;
        k=k+1; 
        if k > maxiter
            info('k>maxiter');
            break
        end
    end    
end