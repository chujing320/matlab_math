function g=gP153(x0)
    x=sym ('x',[8,1]);
    g=zeros(8,1);
    f1=x(1)+x(2);
    f2=x(3)+x(4);
    f3=x(5)*x(1)+x(6)*x(2)-x(7)*x(3)-x(8)*x(4);
    f4=x(7)*x(1)+x(8)*x(2)+x(5)*x(3)+x(6)*x(4);
    f5=x(1)*(x(5)^2-x(7)^2)-2*x(5)*x(3)*x(7)+x(2)*(x(6)^2-x(8)^2)-2*x(8)*x(4)*x(6);
    f6=x(3)*(x(5)^2-x(7)^2)+2*x(5)*x(1)*x(7)+x(4)*(x(6)^2-x(8)^2)+2*x(8)*x(2)*x(6);
    f7=x(1)*x(5)*(x(5)^2-3*x(7)^2)+x(3)*x(7)*(x(7)^2-3*x(5)^2)+x(2)*x(6)*(x(6)^2-3*x(8)^2)+x(4)*x(8)*(x(8)^2-3*x(6)^2);
    f8=x(3)*x(5)*(x(5)^2-3*x(7)^2)-x(1)*x(7)*(x(7)^2-3*x(5)^2)+x(4)*x(6)*(x(6)^2-3*x(8)^2)-x(2)*x(8)*(x(8)^2-3*x(6)^2);
    fun=f1^2+f2^2+f3^2+f4^2+f5^2+f6^2+f7^2+f8^2;
    for i=1:8
        g(i)=subs(diff(fun,x(i)),x,x0); 
    end
end