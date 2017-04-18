function[f0,g0,G0] = Penalty(x0,returnnum)
    if  nargin==1
        returnnum = 3;
    end
    f0 = getPenalty(x0);
    g0= gPenalty(x0);
    %若指定返回的参数个数为2，则无需计算G矩阵
    if returnnum ==2 
        G0 = 0;
        return ;
    end
    G0 = ggPenalty(x0);
end
