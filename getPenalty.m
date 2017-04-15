function [funstr] = getPenalty(n,gamma)
    if nargin==1
        gamma = 10^-5;
    end
    x = sym('x',[n,1]);
    funstr = gamma*sum((x-1).*(x-1))+(sum(x.*x)-1/4)^2;
    
 %   funstr = strcat(num2str(gamma),'*sum((x-1)^2)+(sum(x.*x)-1/4)^2');
end