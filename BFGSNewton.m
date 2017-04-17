function [data_f,data_g, x0 ,feva ] = BFGSNewton(ObjFun,x0,tol,maxiter)
%
% BFGSNewton Newton's Method
%
% Input:
%   ObjFun - input funtion ,choice :Penalty | Chebyquad | p153
%   x0 - x in P, n-vector
%   tol: �����������ֵ
%   maxiter - maximum number of iterations
%
% Output:
% data_f:�洢ÿһ��������fֵ
% daga_g:�洢ÿһ��������gֵ
% x0 :������Ž��xֵ
% feva:�������ô���
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
    
    Rule.crtr='bostwlf';
    Rule.mthd = 'bointrplt33';
    Rule.opt=[1 20 10 0.95 0.05];
    
    k=1;
    x1=x0;
    [n,t] = size(x0);
    [f0 g0]=feval(ObjFun, x0, 2);
    feva = 2;
    H0 = ones(n,n);
    H1 = H0;        
    while norm(g0)>=tol
        %data ��������м�����
        data_f(:,k) = f0;
        data_g (:,k) = g0;
        %��ʼ����
        d = -H0*g0;
        [alaph,info1] = bolinesearch(ObjFun, x0, d, Rule);
        x1 = x0+alaph*d;%������׼��
        [f0 g0]=feval(ObjFun, x1, 2);%���뷵��ֵ����2
        feva = feva+3;
        %������ʽ
        s = x1-x0;
        y = g1-g0;
        H1 = H0+(1+(y'*H0*y)/(y'*s))*((s*s')/(y'*s))-((s*y'*H0+H0*y*s')/(y'*s));
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