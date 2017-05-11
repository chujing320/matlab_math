function [alpha, feva] = mybostwolf(ObjFun, x0, dk, feva ,f0, g0, beta, rho)

if nargin==6 %设置默认值
    beta = 0.55; rho = 0.4;
end 
m = 0; mk = 0;
while m<30 %amijo准则
 xk = x0+beta^m*dk;
 f1 = feval(ObjFun,xk,1);
    if (f1<f0+rho*beta^m*g0'*dk) %&& (abs(g0'*dk) <= -rho*g0'*dk) %armijo & wolf
        mk = m;
        break;
    end
    m = m+1;
    feva = feva+1;
end
% if mk==0
%     mk =1;
% end
mk;
alpha = beta^mk; %如果找不到合适的步长 则mk=0，即alpha=1
end