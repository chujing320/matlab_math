function [ x, ex ] =  DampedNewton( ObjFun,x0,tol,maxiter,varargin)
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
    if strcmp(ObjFun,'Penalty')  
        gamma = 10^-5;
        func = getPenalty(n,gamma);    
    elseif strcmp(ObjFun,'Chebyquad')
         func = getChebyquad(gamma);
    elseif  strcmp(ObjFun,'p153')
        func = getP153(gamma);
    else
         error('DampedNewton: invalid input ObjFun');
    end
    
    k=1;
    x1 = x0;  
    g = jacobian(func);
    g0 = (eval(subs(g,x,x1)))'
    G = hessian(func);
    while norm(g0)>=tol
        f = eval(subs(func,x,x1));
        G0 = eval(subs(G,x,x1))
        d = -g0\G0;
        alaph =1;
        x1 = x1+alaph*d';%ÏßËÑË÷×¼Ôò
        k=k+1;
        if k >maxiter
            info('k>MaxxIter');
            break
        end
        g0 = (eval(subs(g,x,x1)))';
    end            


end


