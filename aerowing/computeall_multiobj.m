function [f1,f2,c1,ceq1,OutputVal] = computeall_multiobj(x)
%% Parameters/Variables assignment
v = x(end);
xDesign = x(1:end-1);

%% Obtain the initial pressures
para.numNodes = 20; % number of points on the airfoil
para.Airfoil = 'NACA'; % it means a NACA airfoil is used
para.NACA = '0012'; % or 0012
para.width = 0.01; % span, unit: m
para.rho = 1.206; % air density
para.U = v;
para.alpha = 1;
para.Re = 2.78e5; % Reynold number
para.Mach = para.U/343;

% XFOIL simulation
FluidOutput = Aerodynamics_XFOIL(para);
Cp = FluidOutput.Cp; % positive: out of the airfoil
P_initial = FluidOutput.P;

% Write pressures in the text file
% update pressures
filename_pressure = 'Pressures.txt';
fileID_pressure = fopen(filename_pressure,'w');
fprintf(fileID_pressure,'%.4f\n',P_initial);
fclose(fileID_pressure);

%% Define the filename for eigenvalue info
filename_output = 'WingOutput.txt';
% Delete if the previous one exists
if exist(filename_output, 'file')==2
    delete(filename_output);
end

%% Read the design variables
filename_input = 'DesignVariables.txt';
% Delete if the previous one exists
if exist(filename_input, 'file')==2
    delete(filename_input);
end
fileDV = fopen(filename_input,'w');
fprintf(fileDV,'%.4f\n',xDesign);
fclose(fileDV);

%% Call out python code and run FEA
system('cd C:\Users\yktsai0121\Desktop\AERO604_abaqus\CaseStudy_AeroWing');
disp('Starting FEA analysis in Abaqus...')
% tic
system('abaqus cae noGUI=wing_assembly.py');
disp('Done!')
% toc

%% Read the data
fileID = fopen(filename_output,'r');
OutputVal_tmp = fscanf(fileID,'%f');
OutputVal = OutputVal_tmp';
fclose(fileID);
fprintf('Design variables:\n')
fprintf('Input: %.4f\n',xDesign)
fprintf('Output: %.4f\n',OutputVal)

if isempty(OutputVal)
    maxStress = NaN;
    Mass = NaN;
    TipDisp = NaN;
    EigenVal = NaN;
    Buckle = 1;
    Yield = 1;
    Design = 1;
else
    maxStress = OutputVal(:,1);
    Mass = OutputVal(:,2);
    TipDisp = OutputVal(:,3);
    EigenVal = OutputVal(:,4);
    Buckle = OutputVal(:,5);
    Yield = OutputVal(:,6);
    Design = OutputVal(:,7);
end

%% Calculate approx cost
ApproxCost = Mass.*(xDesign(end)*417.4851+(1-xDesign(end))*42.58488);
OutputVal = [OutputVal,ApproxCost];

%% Calculate safety factor
safetyFactor = (xDesign(end)*828000000+(1-xDesign(end))*290000000)./maxStress;
OutputVal = [OutputVal,safetyFactor];

%% Save the data to CSV
writematrix([xDesign,v,OutputVal],'data_opt.csv','WriteMode','append')

%% Optimization functions assignment
f1 = Mass; % minimize mass
f2 = -safetyFactor;
% c1(1) = TipDisp-0.025; % TipDisp <= 0.025 m
c1(1) = TipDisp-0.01; % TipDisp <= 0.025 m
% c1(2) = -(safetyFactor-1); % safety factor >= 1
c1(2) = ApproxCost-1e4;
ceq1(1) = Buckle; % no buckle
ceq1(2) = Yield; % no yield
end