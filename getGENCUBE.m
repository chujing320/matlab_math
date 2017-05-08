function fk = getGENCUBE(xk)
    n = length(xk);
    fk = 0;
    for i = 2:n
        fk = fk+(xk(i)-xk(i-1)^3)^2;
    end
    fk = fk + (xk(1) -  1)^2;
end
    


