function output = Aerodynamics_XFOIL(para)
% Input the airfoil coordinates instead of the NACA code.

%% Parameters
width = para.width; % unit: m (not airfoil's, but the wing's)
rho = para.rho; % air density
U = para.U; % speed, unit: m/s
q = 1/2*rho*U^2; % dynamic pressure

%% Xfoil
alpha = para.alpha; % AoA, unit: deg
Re = para.Re; % Reynold number
Mach = para.Mach;
numNodes = para.numNodes;

NACA = para.NACA;
[CL,CD,Cp,X0,Y0] = xfoil(NACA,alpha,Re,Mach,numNodes);

% output.CL_over_CD = CL/CD;
output.CL = CL;
output.CD = CD;
output.RotationMatrix = [cosd(-alpha),-sind(-alpha);sind(-alpha),cosd(-alpha)];
output.X0 = X0;
output.Y0 = Y0;
Rotated_coordinate = output.RotationMatrix*[X0,Y0]';
output.X = Rotated_coordinate(1,:)';
output.Y = Rotated_coordinate(2,:)';

%% Airfoil coordinate
X_mid = (output.X(2:end)+output.X(1:end-1))/2;
Y_mid = (output.Y(2:end)+output.Y(1:end-1))/2;
X_mid = [X_mid;(output.X(1)+output.X(end))/2];
Y_mid = [Y_mid;(output.Y(1)+output.Y(end))/2];
ds = sqrt((X_mid(1:end-1)-X_mid(2:end)).^2+((Y_mid(1:end-1)-Y_mid(2:end)).^2));
ds = [ds;sqrt((X_mid(1)-X_mid(end))^2+((Y_mid(1)-Y_mid(end))^2))];

X_diff = X_mid(2:end)-X_mid(1:end-1);
Y_diff = Y_mid(2:end)-Y_mid(1:end-1);
X_diff = [X_diff;X_mid(1)-X_mid(end)];
Y_diff = [Y_diff;Y_mid(1)-Y_mid(end)];

% slope
slope = atan2(Y_diff,X_diff);
normal_dir = (slope+3*pi/2);

output.normal_dir = normal_dir;

output.dA = ds*width;

output.Cp = Cp;
output.P = q*Cp;
output.L = q*CL;
