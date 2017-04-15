function [funstr] = getPenalty(gamma)
    if nargin==0
        gamma = 10^-5;
    end
    funstr = strcat(num2str(gamma),'*sum((x-1)^2)+(sum(x.*x)-1/4)^2');
end