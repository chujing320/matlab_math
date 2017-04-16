function [data_f,data_g, k,feva] = SR1Newton(ObjFun,x0,tol,maxiter)

    if nargin==2
        tol=1e-8;
        maxiter = 200;
    elseif nargin==3
        maxiter = 200;
    elseif nargin<2 || nargin>4
        err('error input');
    end
    feva =0;
    [n,t] = size(x0);
    x = sym('x',[n,1]);
    k=1;
    x1=x0;
    H0 = ones(n,n);
    H1 = H0
    if strcmp(ObjFun,'Penalty')  
        gamma = 10^-5;
        func = getPenalty(n,gamma);    
        %g0 = gPenalty(x0, gamma);
        g = jacobian(func);
        g0 = (eval(subs(g,x,x0)))'
        feva = feva+1;
        while norm(g0)>=tol
            f = eval(subs(func,x,x0));
            feva = feva+1;
            %data 用来存放中间数据
            data_f(:,k) = f;
            data_g (:,k) = g0;
            %开始迭代
            d = -H0*g0; 
            alpha = 1;
            x1 = x0+alpha*d;%线搜索准则
            %g1 = gPenalty(x1,gamma);
            g1 = (eval(subs(g,x,x1)))';
            feva = feva+1;
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
    elseif strcmp(ObjFun,'Chebyquad')
        g0 = gChebyquad(x0);
        feva = feva+1;
        while norm(g0)>=tol
            f = getChebyquad(x0);
            feva = feva+1;
            %data 用来存放中间数据
            data_f(:,k) = f;
            data_g (:,k) = g0;
            %开始迭代
            d = -H0*g0; 
            alpha = 1;
            x1 = x0+alpha*d;%线搜索准则
            %g1 = gPenalty(x1,gamma);
            g1 = gChebyquad(x1);
            feva = feva+1;
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
    elseif  strcmp(ObjFun,'p153')
        g0 = gP153(x0);
        feva = feva+1;
        while norm(g0)>=tol
            f = getP153(x0);
            feva = feva+1;
            %data 用来存放中间数据
            data_f(:,k) = f;
            data_g (:,k) = g0;
            %开始迭代
            d = -H0*g0; 
            alpha = 1;
            x1 = x0+alpha*d;%线搜索准则
            %g1 = gPenalty(x1,gamma);
            g1 = gP153(x1);
            feva = feva+1;
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
    else
         error('DampedNewton: invalid input ObjFun');
    end
          
        



