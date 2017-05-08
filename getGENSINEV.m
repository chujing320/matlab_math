function f = getGENSINEV(x0)
    c1 = 10^(-4); c2 = 4; f =0;
    n = length(x0);
    for i=1:n-1
        f = (x(i+1)-sin(x(x)))^2/c1 + x(i)^2/c2;
    end
end