function [data_f,data_g, x0, feva ] =  DampedNewton( ObjFun,x0,tol,maxiter, varargin)
%
% DampedNewton Newton's Method
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
        maxiter = 5000;
    elseif nargin==3
        maxiter = 5000;
    elseif nargin<2 || nargin>4
        err('error input');
    end
    
    Rule.crtr='bostwlf';
    Rule.mthd = 'bointrplt33';
    Rule.opt=[1 20 10 0.95 0.05];

    [n,t] = size(x0);
    k=1;
    [f0 g0 G0]=feval(ObjFun, x0, varargin{:});
    feva =3;
    while norm(g0)>=tol
        data_f(:,k) = f0;
        data_g (:,k) = g0;
        %d = -G0^-1 * g0;
        d = -G0\(g0+10^-20);
        [alaph,info1] = bolinesearch(ObjFun, x0, d, Rule);
        if info1(1) == 1%��û���ҵ�����׼��Ĳ���������Ĭ�ϲ���Ϊ1��ţ�ٷ�
           alaph = 1;
        end
        feva = feva+1;
        x0 = x0+alaph*d;%������׼��
        k=k+1;
        if k >maxiter
            info('k>MaxxIter');
            break
        end
        [f0,g0,G0]=feval(ObjFun, x0);
        feva = feva+3;
     end         

end


