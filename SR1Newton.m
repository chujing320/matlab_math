function [x1,k,data] = SR1Newton(ObjFun,x0,tol,maxiter)

    if nargin==2
        tol=1e-8;
        maxiter = 200;
    elseif nargin==3
        maxiter = 200;
    elseif nargin<2 || nargin>4
        err('error input');
   
    [n,t] = size(x0);
    if strcmp(ObjFun,'Penalty')  
        gamma = 10^-5;
        func = getPenalty(gamma);
    elseif strcmp(ObjFun,'Chebyquad')
         func = getChebyquad(gamma);
    elseif  strcmp(ObjFun,'p153')
        func = getP153(gamma);
    else
         error('DampedNewton: invalid input ObjFun');
    end
    k=0;
    x1=x0;
    H0 = ones(n,n);
    H1 = H0;
    %g0 = gPenalty(x0, gamma);
    g0 = jacobian(func);
    while norm(g0)>=tol
        f = eval(subs(func,x,x0));
        %�˴����k/x/f/g ȷ��һ���ÿ��ĸ�ʽ
        perStepPrinf(k,x0,f,g0);
        %data ��������м�����
        data(:,k)=x0;
        %��ʼ����
        d = -H0*g0; 
        x1 = x0+alpha*d;%������׼��
        %g1 = gPenalty(x1,gamma);
        g1 = jacobian(func);
        %������ʽ
        s = x1-x0;
        y = g1-g0;
        H1 = (H0+(s-H0*y)*(s-H0*y)')/(s-H0*y)'*y;
        x0 = x1;
        g0 = g1;
        H0 = H1;
        k=k+1; 
        if k > maxiter
            info('k>maxiter');
            break
        end
    end            
        



