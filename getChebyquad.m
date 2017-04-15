function func=getChebyquad(x_k)
    func=0;
    [n,t] = size(x_k);
    for i = 1:n
        a_i = ChebyshevPoly(i);
        r_i = sum(polyval(a_i,x_k))/n - integral(@(x)polyval(a_i,x),0,1);
        func = func +r_i.^2;
    end
end