function f = parameterfunA(x,p)
% m: number of objective functions
% n: number of design variables
% p: number of parameters
% m = 1; n = length(x); p = 1;
% for j = 1:p
%     f(j) = 1/2*sin(20*x(m+j));
% end
f = zeros(1,p);
for i = 1:p
    f(i) = x(i);
%     f(i) = 1/2*sin(20*x(i))+0.5;
end
% f = 1/2*sin(20*x(1));