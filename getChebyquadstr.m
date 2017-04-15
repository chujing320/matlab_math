function [funstr] = getChebyquadstr(m,n)
%求切比雪夫目标函数 
    if nargin==1
        n = m;
    end

    j=1;
    fstr = '';
    cheby1 = '';
    cheby2 = '';
    for i=1:m
        cheby2 = getChebyploy(i); 
        while j<=n
            cheby1 = getChebyploy(i,j);
            fstr = strcat(fstr,cheby1,'/n-int(',cheby2,'0,1)');
            j = j+1;
        end
    end
    funstr = fstr;
