function [data_f, data_g, x0, feva] = CG2(ObjFun, x0, tol, maxiter)

    if nargin==2
        tol=1e-8;
        maxiter = 2000;
    elseif nargin==3
        maxiter = 2000;
    elseif nargin<2 || nargin>4
        err('error input');
    end
    
    rho = 0.55; sigma = 0.4; k = 1;
    [n,t] = size(x0);
    g0 = ones(n,t); % ��ʼ��g0Ϊ1����
    dk = zeros(n,t); % ��ʼ��dkΪ0����
    [f1,g1] = feval(ObjFun, x0, 2);
    while norm(g1)>=tol
        feva = feva+2;
        data_f(:,k) = f1;
        daga_g(:,k) = g1;
        a = g1'*g1/(g0'*g0); %FR����ϵ��
        b = g1'*(g1-g0)/(g1'*g1); %PRP����ϵ��
        if mod(k-1,50)==0 %�ӵ�0�ο�ʼ�� n�����¿�ʼ׼��
            dk = -g1;
        else
            dk = -g1+max(-a,b)*dk;
        end
        
        f0 = f1; g0 = g1; m = 0; mk = 0;
        while m<20 %amijo׼��
            f1 = feval(ObjFun,x0+rho^m*g1'*dk,1);
            if f1<f0+sigma*rho^m*g1'*dk && 
                mk = m;
                break;
            end
            m = m+1;
            feva = feva+1;
        end
        x0 = x0+rho^mk*dk;
        k = k+1;
        if k > maxiter
            info('k>maxiter');
            break
        end
        [f1,g1] = feval(ObjFun, x0, 2);
    end
end