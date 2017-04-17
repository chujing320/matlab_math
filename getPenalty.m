function [func] = getPenalty(n,gamma)
    if nargin==1
        gamma = 10^-5;
    end
    x = sym('x',[n,1]);
    func = gamma*sum((x-1).*(x-1))+(sum(x.*x)-1/4)^2;
 end