function [f0,g0,G0 ] = S303(x0,returnnum)
   if  nargin==1
        returnnum = 3;
   end

   [f0,g0] = getS303(x0,returnnum);
%    f0 = getS303(x0);
%    g0 = gS303(x0);
   if returnnum == 2 
        G0=0;
        return ;
   end
 %   G0=ggS303(x0);
end
