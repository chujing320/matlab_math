function g = gGENCUBE(x)

[n,t] = size(x);
if n==0 && t==0
    error('error input xk');
end
g = zeros(n,t);
g(1) = 2*x(1)-2-600*x(1)^2*(x(2)-x(1)^3);
for i=2:n-1
    g(i) = 200*(x(i)-x(i-1)^3)-600*x(i)^2*(x(i+1)-x(i)^3);
end
g(n) = 200*(x(n)-x(n-1)^3);
end
    