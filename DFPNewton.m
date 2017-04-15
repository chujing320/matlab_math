function [x1,k,data] = DFPNewton(ObjFun,x0,tol,maxiter)

    if nargin==1
        tol=1e-8;
    end
    
    k=0;
    x1=x0;
    [n.t] = size(x0);
    H0 = ones(n,n);
    H1 = H0;
    if strcmp(ObjFun,'Penalty')
        gamma = 10^-5;
        func = @(x,gamma)gamma*sum((x-1)^2)+(sum(x.*x)-1/4)^2;
        g0 = gPenalty(x0, gamma);
        while norm(g0)>=tol
            f = func(x0,gamma);
            %�˴����k/x/f/g ȷ��һ���ÿ��ĸ�ʽ
            perStepPrinf(k,x0,f,g0);
            %data ��������м�����
            data(:,k)=x0;
            %��ʼ����
            d = -H0*g0;
            x1 = x0+alaph*d;%������׼��
            g1 = gPenalty(x1,gamma);
            %������ʽ
            s = x1-x0;
            y = g1-g0;
            H1 = H0+(s*s')/(s'*y)-(H0*y*y'*H0)/(y'*H0*y);
            x0 = x1;
            g0 = g1;
            H0 = H1;
            k=k+1; 
            if k > maxiter
                info('k>maxiter');
                break
            end
        end            
    elseif strcmp(ObjFun,'Chebyquad')
        func = @(x,i)  
        
    elseif strcmp(ObjFun,'p153')
        
    else
         error('DampedNewton: invalid input ObjFun');
    end
