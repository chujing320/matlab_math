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
    %��ָ�����صĲ�������Ϊ2�����������G����
    if returnnum ==2 
        G = 0;
        return ;
    end
    G = hessian(func);
    G0 = eval(subs(G,x,x0));
end
