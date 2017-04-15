function [funstr] = getChebyploy(n,j)
    %求切比雪夫公式表达
    a = ChebyshevPoly(n); %求切比雪夫的系数
    i=1;
    fstr = '';
    if nargin==1
        while i<=n
        f = strcat('(',num2str(a(i)),')*x^',num2str(n-i+1),'+');
        fstr = strcat(fstr,f);
        i=i+1;
        end
        funstr = strcat(fstr,'(',num2str(a(n+1)),')');
    else
        while i<=n
            f = strcat('(',num2str(a(i)),')*x(',num2str(j),')^',num2str(n-i+1),'+');
            fstr = strcat(fstr,f);
            i=i+1;
        end
        funstr = strcat(fstr,'(',num2str(a(n+1)),')');
    end
end