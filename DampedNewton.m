function [ x, ex ] =  DampedNewton( ObjFun,x0,Step,maxiter,RuleMin,varargin)
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
%   x - aproximation to root
%   ex - error estimate
%
% Example:
%	[ x, ex ] = DampedNewton( 'Penalty', 0, 0.5*10^-5, 10 )
%
% Version:  2017.4.10
% Create:   2017.4.10
% Coder:    Chujing Tan

    if isempty(Step)
        Step = zeros(size(x0));
    end
    [n,t] = size(x0);
    k=0;
    syms x;
    x = x0;
    if strcmp(ObjFun,'Penalty')
        gamma = 10^-5;
        func = @(x,gamma) gamma*sum((x-1)^2)+(sum(x.*x)-1/4)^2;
        g = gPenalty(x, gamma);
        while norm(g)>=RuleMin
            f = func(x,gamma);
            G = ggPenalty(x,gamma);
            d = -g\G;
            x = x+alaph*d;%线搜索准则
            %此处输出k/x/f/g 确定一个好看的格式
            perStepPrinf(k,x,f,g);
            k=k+1;
            if k >maxiter
                info('k>MaxxIter');
                break
            end
            g = gPenalty(x, gamma);
        end            
    elseif strcmp(ObjFun,'Chebyquad')
        j=1; %x的下标
        f =0;
        while j<=n         
            f = f + ChebyshevX(j,x,
        end
        
        
    elseif strcmp(ObjFun,'p153')
        
    else
         error('DampedNewton: invalid input ObjFun');
    end

end


