function [x1,k,data_x0,data_fg] = SR1Newton(ObjFun,x0,tol,maxiter)

    if nargin==2
        tol=1e-8;
        maxiter = 200;
    elseif nargin==3
        maxiter = 200;
    elseif nargin<2 || nargin>4
        err('error input');
    end

    fid = fopen(strcat(ObjFun,'_SR1.txt'),'wt');
    fprintf(fid,'k      f     g\n');
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
    x1=x0;
    H0 = ones(n,n);
    H1 = H0
    %g0 = gPenalty(x0, gamma);
    g = jacobian(func);
    g0 = (eval(subs(g,x,x0)))'
    while norm(g0)>=tol
        f = eval(subs(func,x,x0));
        %�˴����k/x/f/g ȷ��һ���ÿ��ĸ�ʽ
        fprintf(fid,'%8i%8f%8f\n',k,f,g0);
        %data ��������м�����
        data_x0(:,k)=x0;
       % data_fg(:,k)=[f,g0];
        %��ʼ����
        d = -H0*g0; 
        alpha = 1;
        x1 = x0+alpha*d;%������׼��
        %g1 = gPenalty(x1,gamma);
        g1 = (eval(subs(g,x,x1)))';
        %������ʽ
        s = x1-x0
        y = g1-g0
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
        



