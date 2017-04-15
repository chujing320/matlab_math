function [funstr] = getChebyploy(n,j)
    %���б�ѩ��ʽ���
    a = ChebyshevPoly(n); %���б�ѩ���ϵ��
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