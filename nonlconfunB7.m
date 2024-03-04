function [c,ceq] = nonlconfunB7(x,m,p)
n = length(x);
% theta = parameterfunA(x,p);
cm = 0.8;
if p>=1 % parametric
    for i = 1:p
        cm = cm+0.1/p*sin(2*pi*x(i));
    end
    c(1) = -(x(p+m)-cm);
end
% if m>1 % multi-objective
%     f = objectivefunA(x,m,p);
%     for i = 1:(m-1)
%         cm = cm+0.1/(n-1)*sin(2*pi*f(i));
%     end
% end

if m>1 % multi-objective
    cm = 0.85; cm2 = 0.65;
%     for j = 1:p
%         cm = cm+0.1/(n-1)*sin(5*pi*x(j));
%         cm2 = cm+0.1/(n-1)*sin(5*pi*x(j));
%     end
    if p>=1 % parametric
        for i = 1:p
            cm = cm+0.15/(n-1)*sin(10*pi*x(i));
            cm2 = cm2-0.15/(n-1)*sin(10*pi*x(i));
        end
    end
    f = objectivefunB(x,m,p);
    for i = 1:(m-1)
        cm = cm+0.05/(n-1)*sin(10*pi*f(i));
        cm2 = cm2-0.05/(n-1)*sin(10*pi*f(i));
    end
    if p>=1
        c(2) = -cm+x(p+m);
        c(3) = cm2-x(p+m);
    else
        c(1) = -cm+x(p+m);
        c(2) = cm2-x(p+m);
    end
end


ceq = [];
end