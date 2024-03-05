function c = geo_constr(x,h0,L)
hL = x(1);
t = x(2);
l1 = x(3);
l2 = x(4);
r1 = x(5);
r2 = x(6);
tol = 0.0001; % unit: m

upper_line_slope = (hL-h0)/L;
center_line_slope = upper_line_slope/2;

y1 = center_line_slope*l1+h0/2; % y-component of center of first circle
y2 = center_line_slope*l2+h0/2; % y-component of center of second circle

c(1) = (h0-hL)/2/L*l1-h0/2+r1+tol;
c(2) = (h0-hL)/2/L*l2-h0/2+r2+tol;

