function [c,ceq] = nonlconfunA8_new(x,m,p)
n = length(x);
cm_old = -0.1*x(1)+0.6; cm_old2 = 0.4;
if p>=1 % parametric
    for i = 1:p
        cm_old = cm_old+0.3/p*sin(1.5*pi*x(i));
%         cm_old2 = cm_old2+0.25/(n-1)*cos( 2*pi*(x(i)-0.5) );
    end
    if m>1
        f = objectivefunA(x,m,p);
        for j = 1:(m-1)
            cm_old2 = cm_old2+0.2*sin(pi*x(1))/(m-1)*cos( pi/2*f(j) );
%             cm_old2 = cm_old2+0.25/(n-1)*cos( 2*pi*(f(j)-0.5) );
        end
    end
    c(1) = -(x(p+m)-cm_old);
    c(2) = -(x(p+m)-cm_old2);
%     c(2) = -1;
end

if m>1 % multi-objective
%     cm_old = -0.1*x(1)+0.7;
%     if m>1
%         f = objectivefunA(x,m,p);
%         cm_old = cm_old+0.1/(n-1)*sin(1.5*pi*f(i));
%     end
    cm = -0.3*x(1)+0.8; cm2 = -0.3*x(1)+0.7;
    f = objectivefunA(x,m,p);
    if p>=1 % parametric
        for i = 1:p
            cm = cm+0.1/(n-1)*cos(4*pi*x(i));
            cm2 = cm2-0.1/(n-1)*cos(4*pi*x(i));
        end
    end
    
    for i = 1:(m-1)
        cm = cm+0.1/(n-1)*cos(4*pi*f(i));
        cm2 = cm2-0.1/(n-1)*cos(4*pi*f(i));
    end
    if p>=1
        c(3) = -cm+x(p+m);
        c(4) = cm2-x(p+m);
    else
        c(1) = -cm+x(p+m);
        c(2) = cm2-x(p+m);
    end
end
ceq = [];
end