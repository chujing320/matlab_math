function [f,g] = getS303(x0,returnnum)
    [n,t] = size(x0);
    f = 0;
    f1 = 0;
    for i=1:n
        f = f + x(i)^2;
        f1 = f1 + x(i)*i/2;
    end
    f = f + f1^2 +f1^4;
    g = zeros(n,t);
    if returnnum ==1 
        return ;
    end
    for i=1:n
        g(i) = 2*x(i)+i*f1+2*i*f1^3;
    end
end
        
        