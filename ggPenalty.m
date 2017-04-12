function G = ggPenalty(Point,gamma)
%��Penalty���麯���Ķ��׵�G
%PointΪn*1��x����
%GΪn*n�ľ���
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

