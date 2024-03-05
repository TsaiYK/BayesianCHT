function f = mass_fun(x,h0,L,density)

hL = x(1);
t = x(2);
l1 = 0.25; l2 = 0.75;
upper_line_slope = (hL-h0)/L;
h_L_over4 = upper_line_slope*l1+h0;
h_L_3over4 = upper_line_slope*l2+h0;
r1 = x(3)*h_L_over4/2;
r2 = x(4)*h_L_3over4/2;

A1 = pi*r1^2;
A2 = pi*r2^2;
f = t*( (h0+hL)*L/2 - A1 - A2)*density;