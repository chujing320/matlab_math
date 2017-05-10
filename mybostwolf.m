function [alpha, feva] = mybostwolf(ObjFun, x0, dk, feva ,f0, g0, rho, sigma)

if nargin==6 %����Ĭ��ֵ
    rho = 0.45; sigma = 0.5;
end 
m = 0; mk = 0;
while m<30 %amijo׼��
 xk = x0+rho^m*dk;
 f1 = feval(ObjFun,xk,1);
    if (f1<f0+sigma*rho^m*g0'*dk) %&& (abs(g0'*dk) <= -sigma*g0'*dk) %armijo & wolf
        mk = m;
        break;
    end
    m = m+1;
    feva = feva+1;
end
alpha = rho^mk %����Ҳ������ʵĲ��� ��mk=0����alpha=1
end