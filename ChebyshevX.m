function [f] = ChebyshevX(n,x,j)
a = ChebyshevPoly(n); %求切比雪夫的系数
i=1;
f =0;
%构造出切比雪夫多项式 a*x^q
while j<= n
    while i<= n
        f = f + a(i)*x(j)^(n-i+1); 
        i=i+1;
    end
    j=j+1;
end
%构造定积分公式
while 
    
