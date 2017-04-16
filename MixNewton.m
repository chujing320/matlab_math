function [data_f,data_g,k,feva] = MixNewton(ObjFun,x0,tol,maxiter,epslon1,epslon2)

    if nargin==2
        tol=1e-8;
        maxiter = 200;
        epslon1 = 10^-8;
        epslon2 = 10^-8;
    elseif nargin==3
        maxiter = 200;
        epslon1 = 10^-8;
        epslon2 = 10^-8;
    elseif nargin<2
        err('error input');
    end
    
    feva =0;
    [n,t] = size(x0);
    x = sym('x',[n,1]);
    k=1;
    epslon1=1;
    epslon2=1;
    if strcmp(ObjFun,'Penalty')
        gamma = 10^-5;
        func = getPenalty(n,gamma);   
        g = jacobian(func);
        g0 = (eval(subs(g,x,x0)))';
        feva = feva+1;
        G = hessian(func);
        while norm(g0)>=tol
            f = eval(subs(func,x,x0));
            feva = feva+1;
            data_f(:,k) = f;
            data_g (:,k) = g0;
            G0 = eval(subs(G,x,x0));
            feva = feva+1;
            if det(G0)<10^-6 %是奇异矩阵
                d=-g0;                  
            else %G0是非奇异矩阵
                d = -G0^-1 * g0;
                if g0'*d >norm(g0)*norm(d)*epslon1
                    d = -d;
                elseif abs(g0'*d)<=norm(g0)*norm(d)*epslon2
                    d =-g0;
                end
            end    
            %为下一步迭代做准备
            %alaph = linesearch(); %线搜索准则求alaph
            alaph =1;
            x0 = x0+alaph*d; 
            k=k+1;
            if k > maxiter
                info('k>MaxxIter');
                break
            end
           g0 = (eval(subs(g,x,x0)))';
           feva = feva+1;
        end   
    elseif strcmp(ObjFun,'Chebyquad')
          g0 = gChebyquad(x0);
          feva = feva+1;
         while norm(g0)>=tol
             f = getChebyquad(x0);
             feva = feva+1;
             data_f(:,k)=f;
             data_g(:,k)=g0;
             G0=ggChebyquad(x0);
             feva = feva+1;
             if det(G0)<10^-6 %是奇异矩阵
                d=-g0;                  
            else %G0是非奇异矩阵
                d = -G0^-1 * g0;
                if g0'*d >norm(g0)*norm(d)*epslon1
                    d = -d;
                elseif abs(g0'*d)<=norm(g0)*norm(d)*epslon2
                    d =-g0;
                end
             end  
             %为下一步迭代做准备
             %alaph = linesearch(); %线搜索准则求alaph
             alaph =1;
             x0 = x0+alaph*d; 
             k=k+1;
            if k > maxiter
                info('k>MaxxIter');
                break
            end
            g0 = gChebyquad(x0);
            feva = feva+1;
         end
        
    elseif strcmp(ObjFun,'P153')
         g0 = gP153(x0);
         feva = feva+1;
         while norm(g0)>=tol
             f = getP153(x0);
             feva = feva+1;
             data_f(:,k)=f;
             data_g(:,k)=g0;
             G0=ggP153(x0);
             feva = feva+1;
             if det(G0)<10^-6 %是奇异矩阵
                d=-g0;                  
            else %G0是非奇异矩阵
                d = -G0^-1 * g0;
                if g0'*d >norm(g0)*norm(d)*epslon1
                    d = -d;
                elseif abs(g0'*d)<=norm(g0)*norm(d)*epslon2
                    d =-g0;
                end
             end  
             %为下一步迭代做准备
             %alaph = linesearch(); %线搜索准则求alaph
             alaph =1;
             x0 = x0+alaph*d; 
             k=k+1;
            if k > maxiter
                info('k>MaxxIter');
                break
            end
            g0 = gP153(x0);
            feva = feva+1;
         end
    else
         error('DampedNewton: invalid input ObjFun');
    end