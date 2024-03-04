function [c,ceq] = nonlconfunB4(x,m,p)
n = length(x);
cm = 0.75;
if p>=1
    theta = parameterfunA(x,p);
    for i = 1:p
        cm = cm+0.1/p*sin(2*pi*theta(i));%+0.05*sin(6*pi*theta(i));
    end
    c(1) = -(x(p+1)-cm);
else
    f = objectivefunB(x,m,p);
    for i = 1:(m-1)
%         cm = cm+0.1/(n-1)*sin(2*pi*f(i));%+0.05*sin(6*pi*theta(i));
        cm = cm+0.2/(n-1)*sin(pi*(f(i)+0.5));
    end
    c(1) = -(x(m)-cm);
end

% theta = parameterfunA(x,p);
% cm = 0.75;
% for i = 1:p
%     cm = cm+0.1/p*sin(2*pi*theta(i));%+0.05*sin(6*pi*theta(i));
% end
% c(1) = -(x(p+1)-cm);
ceq = [];
end