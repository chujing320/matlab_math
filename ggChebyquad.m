function G=ggChebyquad(x_k)
    G=0;
    [n,t]=size(x_k);
    for i = 1:n
        a_i = ChebyshevPoly(i);
        r_i = sum(polyval(a_i,x_k))/n - integral(@(x)polyval(a_i,x),0,1);
        dT_i = polyval(a_i(1:end-1).*(i:-1:1)',x_k);
        ddT_i = polyval(a_i(1:end-2).*(i:-1:2)'.*(i-1:-1:1)',x_k);
        G = G +dT_i'*dT_i*2/(n^2) + diag(ddT_i*2*r_i/n);
    end
end