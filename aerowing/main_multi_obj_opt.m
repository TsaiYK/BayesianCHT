clear
clc
close all

W_constr_handling = 3;
if W_constr_handling==1
    addpath(genpath('C:\Users\yktsai0121\Documents\GitHub\p3ga_constr_handling\p3ga'))
    rmpath(genpath('C:\Users\yktsai0121\Documents\GitHub\p3ga\p3ga'))
    addpath(genpath('C:\Users\yktsai0121\Documents\GitHub\p3ga_constr_handling\dd_tools'))
    rmpath(genpath('C:\Users\yktsai0121\Documents\GitHub\p3ga\dd_tools'))
    addpath(genpath('C:\Users\yktsai0121\OneDrive - Texas A&M University\2022Fall\IEEE_Evo_Comp_2023\matlab_codes'))
elseif W_constr_handling==2
    addpath(genpath('C:\Users\yktsai0121\Documents\GitHub\p3ga\p3ga'))
    rmpath(genpath('C:\Users\yktsai0121\Documents\GitHub\p3ga_constr_handling\p3ga'))
    addpath(genpath('C:\Users\yktsai0121\Documents\GitHub\p3ga\dd_tools'))
    rmpath(genpath('C:\Users\yktsai0121\Documents\GitHub\p3ga_constr_handling\dd_tools'))
    rmpath(genpath('C:\Users\yktsai0121\OneDrive - Texas A&M University\2022Fall\IEEE_Evo_Comp_2023\matlab_codes'))
elseif W_constr_handling==3
    addpath(genpath('C:\Users\yktsai0121\Documents\GitHub\p3ga\p3ga'))
    rmpath(genpath('C:\Users\yktsai0121\Documents\GitHub\p3ga_constr_handling\p3ga'))
    addpath(genpath('C:\Users\yktsai0121\Documents\GitHub\p3ga\dd_tools'))
    rmpath(genpath('C:\Users\yktsai0121\Documents\GitHub\p3ga_constr_handling\dd_tools'))
    addpath(genpath('C:\Users\yktsai0121\Documents\GitHub\PlatEMO\PlatEMO\Algorithms\Multi-objective optimization\C-TAEA'))
    addpath(genpath('C:\Users\yktsai0121\Documents\GitHub\PlatEMO\PlatEMO\Algorithms\Utility functions'))
    addpath(genpath('C:\Users\yktsai0121\Documents\GitHub\PlatEMO\PlatEMO\Problems'))
elseif W_constr_handling==4
    addpath(genpath('C:\Users\yktsai0121\Documents\GitHub\p3ga\p3ga'))
    rmpath(genpath('C:\Users\yktsai0121\Documents\GitHub\p3ga_constr_handling\p3ga'))
    addpath(genpath('C:\Users\yktsai0121\Documents\GitHub\p3ga\dd_tools'))
    rmpath(genpath('C:\Users\yktsai0121\Documents\GitHub\p3ga_constr_handling\dd_tools'))
    rmpath(genpath('C:\Users\yktsai0121\Documents\GitHub\PlatEMO\PlatEMO\Algorithms\Multi-objective optimization\C-TAEA'))
    addpath(genpath('C:\Users\yktsai0121\Documents\GitHub\PlatEMO\PlatEMO\Algorithms\Multi-objective optimization\c-DPEA'))
    addpath(genpath('C:\Users\yktsai0121\Documents\GitHub\PlatEMO\PlatEMO\Algorithms\Utility functions'))
    addpath(genpath('C:\Users\yktsai0121\Documents\GitHub\PlatEMO\PlatEMO\Problems'))
end
addpath(genpath('C:\Users\yktsai0121\Documents\GitHub\randlc'))
addpath(genpath('C:\Users\yktsai0121\Documents\GitHub\randlc\resources\cprnd')) % for some reason this doesn't work sometimes
addpath(genpath('C:\Users\yktsai0121\Documents\GitHub\EPO'));

set(0,'DefaultTextInterpreter','latex'); % change the text interpreter
set(0,'DefaultLegendInterpreter','latex'); % change the legend interpreter
set(0,'DefaultAxesTickLabelInterpreter','latex'); % change the tick interpreter
set(0, 'DefaultLineLineWidth', 1.5);
set(groot,'defaultAxesXGrid','on')
set(groot,'defaultAxesYGrid','on')
set(0,'defaultAxesFontSize',14)

%%
lb = [2, 0.001, 0.0005, 0.0010, 0.01, 1.0, 0, 40];
ub = [10, 0.01, 0.0015, 0.0025, 0.02, 1.3, 1, 65];

x0 = [10	0.009395008179	0.0005	0.00111	0.01	1.293863856  0, 50];
[x_opt,fval,M] = runobjconstr_multiobj(x0,lb,ub,W_constr_handling);

figure;plot3(fval(:,3),fval(:,1),fval(:,2),'o')
xlabel('Parameter: air speed [m/s]')
ylabel('Objective 1: Mass [kg]')
zlabel('Objective 2: -Safety factor')


%%
clear
clc
close all

W_constr_handling = 4;
if W_constr_handling==1
    addpath(genpath('C:\Users\yktsai0121\Documents\GitHub\p3ga_constr_handling\p3ga'))
    rmpath(genpath('C:\Users\yktsai0121\Documents\GitHub\p3ga\p3ga'))
    addpath(genpath('C:\Users\yktsai0121\Documents\GitHub\p3ga_constr_handling\dd_tools'))
    rmpath(genpath('C:\Users\yktsai0121\Documents\GitHub\p3ga\dd_tools'))
    addpath(genpath('C:\Users\yktsai0121\OneDrive - Texas A&M University\2022Fall\IEEE_Evo_Comp_2023\matlab_codes'))
elseif W_constr_handling==2
    addpath(genpath('C:\Users\yktsai0121\Documents\GitHub\p3ga\p3ga'))
    rmpath(genpath('C:\Users\yktsai0121\Documents\GitHub\p3ga_constr_handling\p3ga'))
    addpath(genpath('C:\Users\yktsai0121\Documents\GitHub\p3ga\dd_tools'))
    rmpath(genpath('C:\Users\yktsai0121\Documents\GitHub\p3ga_constr_handling\dd_tools'))
    rmpath(genpath('C:\Users\yktsai0121\OneDrive - Texas A&M University\2022Fall\IEEE_Evo_Comp_2023\matlab_codes'))
elseif W_constr_handling==3
    addpath(genpath('C:\Users\yktsai0121\Documents\GitHub\p3ga\p3ga'))
    rmpath(genpath('C:\Users\yktsai0121\Documents\GitHub\p3ga_constr_handling\p3ga'))
    addpath(genpath('C:\Users\yktsai0121\Documents\GitHub\p3ga\dd_tools'))
    rmpath(genpath('C:\Users\yktsai0121\Documents\GitHub\p3ga_constr_handling\dd_tools'))
    addpath(genpath('C:\Users\yktsai0121\Documents\GitHub\PlatEMO\PlatEMO\Algorithms\Multi-objective optimization\C-TAEA'))
    addpath(genpath('C:\Users\yktsai0121\Documents\GitHub\PlatEMO\PlatEMO\Algorithms\Utility functions'))
    addpath(genpath('C:\Users\yktsai0121\Documents\GitHub\PlatEMO\PlatEMO\Problems'))
elseif W_constr_handling==4
    addpath(genpath('C:\Users\yktsai0121\Documents\GitHub\p3ga\p3ga'))
    rmpath(genpath('C:\Users\yktsai0121\Documents\GitHub\p3ga_constr_handling\p3ga'))
    addpath(genpath('C:\Users\yktsai0121\Documents\GitHub\p3ga\dd_tools'))
    rmpath(genpath('C:\Users\yktsai0121\Documents\GitHub\p3ga_constr_handling\dd_tools'))
    rmpath(genpath('C:\Users\yktsai0121\Documents\GitHub\PlatEMO\PlatEMO\Algorithms\Multi-objective optimization\C-TAEA'))
    addpath(genpath('C:\Users\yktsai0121\Documents\GitHub\PlatEMO\PlatEMO\Algorithms\Multi-objective optimization\c-DPEA'))
    addpath(genpath('C:\Users\yktsai0121\Documents\GitHub\PlatEMO\PlatEMO\Algorithms\Utility functions'))
    addpath(genpath('C:\Users\yktsai0121\Documents\GitHub\PlatEMO\PlatEMO\Problems'))
end
addpath(genpath('C:\Users\yktsai0121\Documents\GitHub\randlc'))
addpath(genpath('C:\Users\yktsai0121\Documents\GitHub\randlc\resources\cprnd')) % for some reason this doesn't work sometimes
addpath(genpath('C:\Users\yktsai0121\Documents\GitHub\EPO'));

set(0,'DefaultTextInterpreter','latex'); % change the text interpreter
set(0,'DefaultLegendInterpreter','latex'); % change the legend interpreter
set(0,'DefaultAxesTickLabelInterpreter','latex'); % change the tick interpreter
set(0, 'DefaultLineLineWidth', 1.5);
set(groot,'defaultAxesXGrid','on')
set(groot,'defaultAxesYGrid','on')
set(0,'defaultAxesFontSize',14)

%%
lb = [2, 0.001, 0.0005, 0.0010, 0.01, 1.0, 0, 40];
ub = [10, 0.01, 0.0015, 0.0025, 0.02, 1.3, 1, 65];

x0 = [10	0.009395008179	0.0005	0.00111	0.01	1.293863856  0, 50];
[x_opt,fval,M] = runobjconstr_multiobj(x0,lb,ub,W_constr_handling);

figure;plot3(fval(:,3),fval(:,1),fval(:,2),'o')
xlabel('Parameter: air speed [m/s]')
ylabel('Objective 1: Mass [kg]')
zlabel('Objective 2: -Safety factor')
