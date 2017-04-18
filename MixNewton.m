function [data_f,data_g,x0,feva] = MixNewton(ObjFun,x0,tol,maxiter,epslon1,epslon2)

% MixNewton Newton's Method
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
%	[ x, ex ] = MixNewton( 'Penalty', [1,1,1,1]' )
%
% Version:  2017.4.10
% Create:   2017.4.10
% Coder:    Chujing Tan

    if nargin==2
        tol=1e-8;
        maxiter =  500;
        epslon1 = 10^-8;
        epslon2 = 10^-8;
    elseif nargin==3
        maxiter = 500;
        epslon1 = 10^-8;
        epslon2 = 10^-8;
    elseif nargin==4
        epslon1 = 10^-8;
        epslon2 = 10^-8;
    elseif nargin<2
        err('error input');
    end
    
    Rule.crtr='bostwlf';
    Rule.mthd = 'bointrplt33';
    Rule.opt=[1 20 10 0.95 0.05];
    [f0 g0 G0]=feval(ObjFun, x0);
    feva =3;
    [n,t] = size(x0);
    k=1;
    f1 =0;
    while norm(g0)>=tol
   % while abs(f0-f1)>=tol
        f0 =f1;
        data_f(:,k) = f0;
        data_g (:,k) = g0;
        if abs(det(G0))<10^(-20) %���������
            d=-g0;                  
        else %G0�Ƿ��������
            d = -G0^-1 * g0+10^-20;
            if g0'*d >norm(g0)*norm(d)*epslon1
                d = -d;
            elseif abs(g0'*d)<=norm(g0)*norm(d)*epslon2
                d =-g0;
            end
        end    
        %Ϊ��һ��������׼��
        [alaph,info1] = bolinesearch(ObjFun, x0, d, Rule);
        if info1(1) == 1%��û���ҵ�����׼��Ĳ���������Ĭ�ϲ���Ϊ1��ţ�ٷ�
            StepSize = 1;
        end
        feva = feva+1;
        x0 = x0+alaph*d; 
        k=k+1;
        if k > maxiter
            info('k>MaxxIter');
            break
        end
        [f1,g0,G0]=feval(ObjFun, x0);
        feva = feva+3;
    end   
end