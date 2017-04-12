function G = ggPenalty(Point,gamma)
%求Penalty检验函数的二阶导G
%Point为n*1的x矩阵
%G为n*n的矩阵
[n,t] = size(Point);
if n==0 && n==0
    error('error input Point');
end
G=zeros(n,n);
Gfunc1 = @(x,i) 8*x(i,1)*x;
Gfunc2 = @(x,i,gamma) 2*gamma+4*(sum(x.*x)+2*x(i,1)*x(i,1))-1;
i=1;
while i<=n
    G(i,:) = Gfunc1(Point,i);
    G(i,i) = Gfunc2(Point,i,gamma);
    i=i+1;
end

