function [c,ceq] = nonlconfunA3(x,m,p)
n = length(x);
cm = 0.6;
if p>=1
    theta = parameterfunA(x,p);
    for i = 1:p
        cm = cm+0.2/p*sin(1.5*pi*theta(i));
    end
    c(1) = -(x(p+1)-cm);
else
    f = objectivefunA(x,m,p);
    for i = 1:(m-1)
        cm = cm+0.2/(n-1)*sin(pi*(f(i)+0.5));
%         cm = cm-0.2*f(i)+0.2*(n-1)*sin(2*pi*f(i));
    end
    c(1) = -(x(m)-cm);
end
ceq = [];
end