function [x,ex] = MixNewton(ObjFun,Point,Step.MaxIter,RuleMin,varargin)

    if isempty(Step)
        Step = zeros(size(Point));
    end
    
    k=0;
    x= Point;
    epslon1=
    epslon2=
    if strcmp(ObjFun,'Penalty')
        gamma = 10^-5;
        func = @(x,gamma)gamma*sum((x-1)^2)+(sum(x.*x)-1/4)^2;
        g = gPenalty(x, gamma);
        while norm(g)>=RuleMin
            f = func(x,gamma);
            G = ggPenalty(x,gamma);
            if det(G)<10^-8 %是奇异矩阵
                d=-g;                  
            else %G是非奇异矩阵
                d = -g\G;
                if g'*d >norm(g)*norm(d)*epslon1
                    d = -d;
                elseif g'*d<=norm(g)*norm(d)*epslon2
                    d =-g;
                end
            end
          
            %此处输出k/x/f/g 确定一个好看的格式
            perStepPrinf(k,x,f,g);
            %为下一步迭代做准备
            alaph = linesearch(); %线搜索准则求alaph
            x = x+alaph*d; 
            k=k+1;
            if k >MaxIter
                info('k>MaxxIter');
                break
            end
            g = gPenalty(x, gamma);
        end   
    elseif strcmp(ObjFun,'Chebyquad')
        
        
    elseif strcmp(ObjFun,'p153')
        
    else
         error('DampedNewton: invalid input ObjFun');
    end