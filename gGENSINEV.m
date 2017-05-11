function g = gGENSINEV(x)
    c1 = 10^(-4); c2 = 4;
    [n,t] = size(x);
    g = zeros(n,t);
    g(1) = -2/c1*(x(2)-sin(x(1)))*cos(x(1)) + 2/c2*x(1);
    for i=2:n-1
        g(i) = 2/c1*(x(i)-sin(x(i-1)))-2/c1*(x(i+1)-sin(x(i)))*cos(x(i)) + 2/c2*x(i);
    end
    g(n) = 2/c1*(x(n)-sin(x(n-1)));
end