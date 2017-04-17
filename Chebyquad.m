function[f0,g0,G0] = Chebyquad(x0,returnnum)
    if  nargin==1
        returnnum = 3;
    end
    f0 = getChebyquad(x0);
    g0 = gChebyquad(x0);
  if returnnum ==2 
        G0=0;
        return ;
    end
    G0=ggChebyquad(x0);
end
