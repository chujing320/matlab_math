function [f0,g0,G0] = Penalty(x0,returnnum)
    if  nargin==1
        returnnum = 3;
    end
    n = length(x0);   
    x = sym('x',[n,1]);
    gamma  = 10^-5; 
    func = getPenalty(n,gamma); 
    f0 = eval(subs(func,x,x0)); 
    g = jacobian(func);
    g0 = (eval(subs(g,x,x0)))';
    %若指定返回的参数个数为2，则无需计算G矩阵
    if returnnum ==2 
        G = 0;
        return ;
    end
    G = hessian(func);
    G0 = eval(subs(G,x,x0));
end
