function [funstr] = getChebyquad(m,n)
%求切比雪夫目标函数
j=1;
fstr = '';
while j<=m
    cheby1 = getChebyploy(i,j);
    cheby2 = getChebyploy(i);
    fstr = strcat(fstr,cheaby1,'/n-int(',cheby2,'0,1)');
    j = j+1;
end
