function [f0, g0, G0] = GENCUBE(x0,returnnum)
    if  nargin==1
        returnnum = 3;
    end
    f0 = getGENCUBE(x0);
    if returnnum == 1;
        g0 = 0; G0 = 0;
        return ;
    end
    g0 = gGENCUBE(x0);
    if returnnum == 2
        G0 = 0;
        return;
    end
  %  G0 = ggGENCUBE(x0);
end


