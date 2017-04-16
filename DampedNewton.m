function [data_f,data_g, k, feva ] =  DampedNewton( ObjFun,x0,tol,maxiter,varargin)
%
% NEWTON Newton's Method
%   Newton's method for finding successively better approximations to the 
%   zeroes of a real-valued function.
%
% Input:
%   ObjFun - input funtion ,choice :Penalty | Chebyquad | p153
%   Point - x in P, n-vector
%   Step:     d in P, n-vector
%   nmax - maximum number of iterations
%
% Output:
%   x - aproximation to rootz
%   ex - error estimate
%
% Example:
%	[ x, ex ] = DampedNewton( 'Penalty', 0, 0.5*10^-5, 10 )
%
% Version:  2017.4.10
% Create:   2017.4.10
% Coder:    Chujing Tan
    if nargin==2
        tol=1e-8;
        maxiter = 200;
    elseif nargin==3
        maxiter = 200;
    elseif nargin<2 || nargin>4
        err('error input');
    end
    
    [n,t] = size(x0);
    x = sym('x',[n,1]);
    k=1;
    feva =0;
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
            d = -G0^-1 * g0;
            alaph =1;
            x0 = x0+alaph*d;%线搜索准则
            k=k+1;
            if k >maxiter
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
             d = -G0^-1 * g0;
             alaph = 1;
             x0 = x0+alaph*d;%线搜索准则
             k=k+1;
            if k >maxiter
                info('k>MaxxIter');
                break
            end
            g0 = gChebyquad(x0);
            feva = feva+1;
         end
             
    elseif  strcmp(ObjFun,'P153')
         g0 = gP153(x0);
         feva = feva+1;
         while norm(g0)>=tol
             f = getP153(x0);
             feva = feva+1;
             data_f(:,k)=f;
             data_g(:,k)=g0;
             G0=ggP153(x0);
             feva = feva+1;
             d = -G0^-1 * g0;
             alaph = 1;
             x0 = x0+alaph*d;%线搜索准则
             k=k+1;
            if k >maxiter
                info('k>MaxxIter');
                break
            end
            g0 = gP153(x0);
            feva = feva+1;
         end
    else
         error('DampedNewton: invalid input ObjFun');
    end
     
   


end


