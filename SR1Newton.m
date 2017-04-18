function [data_f,data_g, x0,feva, g0] = SR1Newton(ObjFun,x0,tol,maxiter)

%SR1Newton Newton's Method
%
% Input:
%   ObjFun - input funtion ,choice :Penalty | Chebyquad | p153
%   x0 - x in P, n-vector
%   tol: 允许的最大误差值
%   maxiter - maximum number of iterations
%
% Output:
% data_f:存储每一步迭代的f值
% daga_g:存储每一步跌倒的g值
% x0 :求得最优解的x值
% feva:函数调用次数
% Example:
%	[ x, ex ] = SR1Newton( 'Penalty', [1,1,1,1]' )
%
% Version:  2017.4.10
% Create:   2017.4.10
% Coder:    Chujing Tan

    if nargin==2
        tol=1e-8;
        maxiter = 10000;
    elseif nargin==3
        maxiter = 10000;
    elseif nargin<2 || nargin>4
        err('error input');
    end
    [n,t] = size(x0);
    
    Rule.crtr='bostwlf';
    Rule.mthd = 'bointrplt33';
    Rule.opt=[1 20 10 0.95 0.05];

    k=1;
    x1=x0;
    H0 = ones(n,n);
    H1 = H0;
    [f0 g0]=feval(ObjFun, x0, 2);
    feva = 2;
    while norm(g0)>=tol
        %data 用来存放中间数据
        data_f(:,k) = f0;
        data_g (:,k) = g0;
        %开始迭代
        d = -H0*g0; 
        [alaph,info1] = bolinesearch(ObjFun, x0, d, Rule);
        if info1(1) == 1%若没有找到满足准则的步长，则用默认步长为1的牛顿法
            StepSize = 1;
        end
        x1 = x0+alaph*d;%线搜索准则
        [f0 g1]=feval(ObjFun, x1, 2);%传入返回值个数2
        feva = feva+2;
        %修正公式
        s = x1-x0;
        y = g1-g0;
        H1 = H0+((s-H0*y)*(s-H0*y)')/((s-H0*y)'*y+10^(-20));
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