function [f] = ChebyshevX(n,x,j)
a = ChebyshevPoly(n); %���б�ѩ���ϵ��
i=1;
f =0;
%������б�ѩ�����ʽ a*x^q
while j<= n
    while i<= n
        f = f + a(i)*x(j)^(n-i+1); 
        i=i+1;
    end
    j=j+1;
end
%���춨���ֹ�ʽ
while 
    
