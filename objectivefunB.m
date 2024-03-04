function f = objectivefunB(x,m,p)
% m: number of objective functions
% n: number of design variables
nvars = length(x);
if m==1
    theta = parameterfunA(x,p);
    g = g_fun(x,m,p,nvars);
    h = 3*(m+p); h_term2 = 0;
    for j = 1:p
        h_term2 = h_term2+2*theta(j)+sin(3*pi*theta(j));
    end
    h = h+h_term2/5/g;
    f = g*h;
else
    theta = parameterfunA(x,p);
    h = 0;
    for j = 1:(m-1)
%         g(j) = g_fun(x,m,p,p+j);
%         h = g(j)*sqrt(1-(f_bar/g(j)/5)^2);
%         f(j) = g(j)*h;
%         f(j) = 1/2*sin(20*x(p+j));
%         f(j) = 2*(x(p+j)-0.5)^2;
        f(j) = x(p+j);
    end
    
    g1 = g_fun(x,m,p,nvars);
    h = 3*(m+p); h_term2 = 0;
    for k = 1:p
        h_term2 = h_term2+2*theta(k)+sin(3*pi*theta(k));
    end
%     for j = 1:(m-1)
%         h_term2 = h_term2-(2*f(j)+sin(3*pi*f(j)));
%     end
    for j = 1:(m-1)
        h_term2 = h_term2+10*sqrt(1-f(j)^2);
    end

    h = h+h_term2/5/g1;
    f(m) = g1*h;
end

function g = g_fun(x,m,p,n)
    x_opt = 0.5;
    for i = 1:(n-1)
        x_opt = x_opt+0.01/(n-1)*sin(5*pi*x(i));
    end
    g = 1;
    g = g+5*( x(n)-x_opt )^2;
%         g = 1+sin(5*pi*x(n));
end

end