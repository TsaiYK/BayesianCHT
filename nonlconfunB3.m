function [c,ceq] = nonlconfunB3(x,m,p)
n = length(x);
% theta = parameterfunA(x,p);
% cm = 0.6;
% for i = 1:p
%     cm = cm+0.3/p*sin(1.5*pi*theta(i));
% end
% c(1) = -(x(p+1)-cm);
% % c(2) = x(p+1)-cm2;
% ceq = [];

cm = 0.6;
if p>=1
    theta = parameterfunA(x,p);
    for i = 1:p
        cm = cm+0.3/p*sin(1.5*pi*theta(i));
    end
    c(1) = -(x(p+1)-cm);
else
    f = objectivefunB(x,m,p);
    for i = 1:(m-1)
        cm = cm+0.3/(n-1)*sin(pi*(f(i)+0.5));
    end
    c(1) = -(x(m)-cm);
end
ceq = [];

end