function [f0] = getPenalty(x0)
    n = length(x0);
    if nargin==1
        gamma = 10^-5;
    end
    x = sym('x',[n,1]);
    func = gamma*sum((x-1).*(x-1))+(sum(x.*x)-1/4)^2;
    f0 = eval(subs(func,x,x0));
 end