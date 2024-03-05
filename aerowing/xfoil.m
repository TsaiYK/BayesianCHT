function [CL,CD,Cp,X_0,Y_0] = xfoil(NACA,AoA,Re,Mach,numNodes)
% CALLING XFOIL FROM MATLAB
% Written by: JoshTheEngineer
% YouTube: www.youtube.com/joshtheengineer
% Website: www.joshtheengineer.com
% Started: 01/01/19
% Updated: 01/01/19 - Started code
%                   - Works as expected

% clear;
% clc;
% tic
% 
% NACA       = '2414';                                                        % NACA 4-digit airfoil [str]
% AoA        = '0';                                                           % Angle of attack [deg]
% numNodes   = '150';                                                          % Panel nodes [#]
% Re = 2.78e5; % Reynold number
% Mach = 10/343;
AoA = num2str(AoA);
numNodes = num2str(numNodes);
saveFlnmAF = 'Save_Airfoil.txt';                                            % Airfoil coordinates filename
saveFlnmCp = 'Save_Cp.txt';                                                 % Pressure coefficient filename
saveFlnmLD = 'Save_LD.pol';                                                 % Lift and Drag coefficients

% Delete files if they exist
if (exist(saveFlnmAF,'file'))
    delete(saveFlnmAF);
end
if (exist(saveFlnmCp,'file'))
    delete(saveFlnmCp);
end
if (exist(saveFlnmLD,'file'))
    delete(saveFlnmLD);
end

% Create the airfoil
fid = fopen('xfoil_input.txt','w');
fprintf(fid,['NACA ' NACA '\n']);
fprintf(fid,'PPAR\n');
fprintf(fid,['N ' numNodes '\n']);
fprintf(fid,'\n\n');

% Save the airfoil data points
fprintf(fid,['PSAV ' saveFlnmAF '\n']);

% Find the Cp vs. X plot
fprintf(fid,'OPER\n');
fprintf(fid,'iter\n');
fprintf(fid,'500\n');
fprintf(fid,['Re ' num2str(Re) '\n']);
fprintf(fid,['Mach ' num2str(Mach) '\n']);
fprintf(fid,'visc \n');
fprintf(fid,['Alfa ' AoA '\n']);
fprintf(fid,['CPWR ' saveFlnmCp '\n']);
fprintf(fid,['PACC ' saveFlnmLD '\n']);
fprintf(fid,[saveFlnmLD '\n']);
fprintf(fid,['\n']);
fprintf(fid,['Alfa ' AoA '\n']);
fprintf(fid,['PACC ' saveFlnmLD '\n']);

% Close file
fclose(fid);

% Run XFoil using input file
cmd = 'xfoil.exe < xfoil_input.txt';
[status,result] = system(cmd);

%% READ DATA FILE: AIRFOIL

% saveFlnmAF = 'Save_Airfoil.txt';                                            % File name
fidAirfoil = fopen(saveFlnmAF);                                             % Open file for reading

dataBuffer = textscan(fidAirfoil,'%f %f','CollectOutput',1,...              % Read data from file
                                 'Delimiter','','HeaderLines',0);
fclose(fidAirfoil);                                                         % Close file
delete(saveFlnmAF);

% Separate boundary points
XB = dataBuffer{1}(:,1);                                                    % Boundary point X-coordinate
YB = dataBuffer{1}(:,2);                                                    % Boundary point Y-coordinate

%% READ DATA FILE: PRESSURE COEFFICIENT

% saveFlnmCp = 'Save_Cp.txt';                                                 % File name
fidCP = fopen(saveFlnmCp);                                                  % Open file for reading
dataCp = textscan(fidCP,'%f %f %f','HeaderLines',3,...                  % Ready data from file
                            'CollectOutput',1,...
                            'Delimiter','');
fclose(fidCP);                                                              % Close file
delete(saveFlnmCp);

% Separate Cp data
X_0  = dataCp{1,1}(:,1);                                                % X-coordinate
Y_0  = dataCp{1,1}(:,2);                                                % Y-coordinate
Cp = dataCp{1,1}(:,3);                                                % Cp data
Cp(1) = 0;
Cp(end) = 0;

%% READ DATA FILE: LIFT & DRAG COEFFICIENT

% saveFlnmLD = 'Save_LD.pol';                                                 % File name
fidLD = fopen(saveFlnmLD);                                                  % Open file for reading
dataLD = textscan(fidLD,'%f','HeaderLines',12);                         % Ready data from file
fclose(fidLD);                                                              % Close file
delete(saveFlnmLD);

% Separate Cp data
CL  = dataLD{1,1}(2);                                                   % Lift
CD  = dataLD{1,1}(4);                                                   % Drag

%% PLOT DATA

% % Split airfoil into (U)pper and (L)ower
% XB_U = XB(YB >= 0);
% XB_L = XB(YB < 0);
% YB_U = YB(YB >= 0);
% YB_L = YB(YB < 0);
% 
% % Split Xfoil results into (U)pper and (L)ower
% Cp_U = Cp(YB >= 0);
% Cp_L = Cp(YB < 0);
% X_U  = X_0(YB >= 0);
% X_L  = X_0(YB < 0);
% 
% % Plot: Airfoil
% figure(1);
% cla; hold on; grid off;
% set(gcf,'Color','White');
% set(gca,'FontSize',12);
% plot(XB_U,YB_U,'b.-');
% plot(XB_L,YB_L,'r.-');
% xlabel('X Coordinate');
% ylabel('Y Coordinate');
% axis equal;
% 
% % Plot: Pressure coefficient
% figure(2);
% cla; hold on; grid on;
% set(gcf,'Color','White');
% set(gca,'FontSize',12);
% plot(X_U,Cp_U,'bo-','LineWidth',2);
% plot(X_L,Cp_L,'ro-','LineWidth',2);
% xlabel('X Coordinate');
% ylabel('Cp');
% ylim('auto');
% set(gca,'Ydir','reverse')
% toc