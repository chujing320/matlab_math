function [f0,g0,G0 ] = GENSINEV(x0,returnnum)
   if  nargin==1
        returnnum = 3;
   end

   f0 = getGENSINEV(x0);
   g0 = gGENSINEV(x0);
  if returnnum == 2 
        G0=0;
        return ;
    end
    G0=ggGENSINEV(x0);
end
