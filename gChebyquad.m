function g=gChebyquad(x_k)
    g=0;
    n=length(x_k);
    for i = 1:n
        a_i = chebyshev(i);
        r_i = sum(polyval(a_i,x_k))/n - integral(@(x)polyval(a_i,x),0,1);
        dT_i = polyval(a_i(1:end-1).*(i:-1:1)',x_k);
        g = g +dT_i.*r_i*2/n;
    end
end