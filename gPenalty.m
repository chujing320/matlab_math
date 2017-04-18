function g = gPenalty(Point)
%求Penalty检验函数的一阶导g
%Point为n*1的x矩阵

[n,t] = size(Point);
if n==0 && t==0
    error('error input Point');
end
gamma = 10^-5;
g = zeros(n,t);
gfunc = @(x,i,gamma) 2*gamma*(x(i,1)-1)+4*(sum(x.*x)-1/4)*x(i,1);
i=1;
while i<=n
    g(i,1) = gfunc(Point,i,gamma);
    i=i+1;
end

