% This script is to run one-time parametric optimization for the numerical
% test problem. The objective function is based on MATLAB tool function
% 'peaks()' defined explicitly by 'myfun', and the constraint is a
% nonlinear function defined by 'nonlcon_fun'. There is only one design
% variable and one parameter in order to easily visualize the results.
% ------------------------------------------------------------------------
% Before running the code, please make sure you include the 'p3ga'-related
% files and add their paths shown in the second section.

clear
clc
close all

addpath(genpath(pwd));

set(0,'DefaultTextInterpreter','latex'); % change the text interpreter
set(0,'DefaultLegendInterpreter','latex'); % change the legend interpreter
set(0,'DefaultAxesTickLabelInterpreter','latex'); % change the tick interpreter
set(0, 'DefaultLineLineWidth', 1);
set(groot,'defaultAxesXGrid','on')
set(groot,'defaultAxesYGrid','on')
set(0,'defaultAxesFontSize',14)

%% Section 1: Choose w/ or w/o constraint handling technique
W_constr_handling = 3;

%% Section 2: Add paths
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
end
addpath(genpath('C:\Users\yktsai0121\Documents\GitHub\randlc'))
addpath(genpath('C:\Users\yktsai0121\Documents\GitHub\randlc\resources\cprnd')) % for some reason this doesn't work sometimes
addpath(genpath('C:\Users\yktsai0121\Documents\GitHub\EPO'));

%% Section 3: Setting for P3GA
h0 = 0.5; L = 1;
density = 2720; % unit: kg/m^3
yielding_stress = 200e6;
Disp_allowed = 4e-3;

obj_fun = @(x) mass_fun(x,h0,L,density);    % objective function
fun1and2 = @(x) [x(end), obj_fun(x)]; % all attributes (including parameter or parameter function)
par = [1];                          % the parameter function is in the first column of fun1and2
dom = [2];                          % the dominator function is in the second column of fun1and2
% Bound constraints on X = [hL, t, fd1, fd2, load]
lb = [0.05,  0.5, 0.1, 0.1,  2.5];
ub = [ 0.5,  2.5, 0.9, 0.9, 10.0];
nvars = length(lb);                 % the number of design variables
A = []; b = []; Aeq = []; beq = []; % linear inequality and equality constraints
Generations = 50;                   % maximum number of generations
PopulationSize = 50;                % maximum populations
ref_bnds = [2.5,10;300,3000];     % reference bounds for HVI measure
nonlcon = @(x) abaqus_constr(x,h0,L,yielding_stress,Disp_allowed); % nonlinear constraint function
classif_err_allowance = 0.2;
TolCon = 1e-6;

%% Section 4: Define 'options'
if W_constr_handling==2
    % Initialize population
    [X_new, feasible, initial, initialfeasible, cEval] = Initialize_LHS_GP(lb,ub,...
        A,b,Aeq,beq,nonlcon,PopulationSize,nvars,'DontPlot',classif_err_allowance,TolCon);
    % Define options
    options = p3gaoptimset('Generations',Generations,'PopulationSize',PopulationSize,...
        'ViewProgress',true,'InitialPopulation',X_new,...
        'CrossoverFraction',0.8,'MutationFraction',0.05);
    options.initial = initial;
    options.initialfeasible = initialfeasible;
    options.GenerationPlots = 1;
    options.classif_err_allowance = classif_err_allowance;
    options.cEval = cEval;
    options.TolCon = TolCon;
else
    % Define options
    options = p3gaoptimset('Generations',Generations,'PopulationSize',PopulationSize,...
        'ViewProgress',true,'CrossoverFraction',0.8,'MutationFraction',0.05);
end

% Save data by generation in a specified directionary
options.Log = true;
options.GenerationData = true;
if options.GenerationData
    if W_constr_handling==2
        options.SaveLoc = strcat('data\test1_w');
    elseif W_constr_handling==3
        options.SaveLoc = strcat('data\test3_w_CTAEA');
    else
        options.SaveLoc = strcat('data\test1_wo');
    end
    if ~exist(options.SaveLoc, 'dir')
        mkdir(options.SaveLoc)
    end
    options.GenerationDataIncrement = 1;
end
options.Dominance = 'hch'; % determine if you want to use HCH-based dominance
options.hvi_options.phi = 45; % Hypercone angle
options.ref_bnds = ref_bnds;

% if W_constr_handling
%     % Save the plot for initial population
%     figure(1);
%     title(['Initial Population'])
%     saveas(figure(1),strcat(options.SaveLoc,'\pop_sample',num2str(0),'.png'))
%     % close(gcf)
% end

%% Section 5: Run P3GA
% Minimize using P3GA
if W_constr_handling==3
    [x_opt,fval,M] = p3ga_TAEA(fun1and2,dom,par,nvars,A,b,Aeq,beq,lb,ub,nonlcon,options);
else
    [x_opt,fval,M] = p3ga(fun1and2,dom,par,nvars,A,b,Aeq,beq,lb,ub,nonlcon,options);
end

figure(gcf); xlabel('Generation','Interpreter','Latex'); ylabel('pHVI','Interpreter','Latex');
saveas(figure(gcf),strcat(options.SaveLoc,'\p3ga_pHVI.fig'))
close(gcf)

figure(gcf); xlabel('$\theta$','Interpreter','Latex'); ylabel('$J$','Interpreter','Latex');
xlim([lb(end),ub(end)])
saveas(figure(gcf),strcat(options.SaveLoc,'\p3ga_result.fig'))
close(gcf)

pHVI_value = M.hvi(2,end);

%% Section 8: Save data
save(strcat(options.SaveLoc,'\matlab.mat'))

if W_constr_handling==2
    figure
    plot(1:Generations,M.exp_classif_err)
    xlabel('Generations'); ylabel('Expected value of classification errors')
end

%% Section 1: Choose w/ or w/o constraint handling technique
clear
clc
close all
W_constr_handling = 4;

%% Section 2: Add paths
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

%% Section 3: Setting for P3GA
h0 = 0.5; L = 1;
density = 2720; % unit: kg/m^3
yielding_stress = 200e6;
Disp_allowed = 4e-3;

obj_fun = @(x) mass_fun(x,h0,L,density);    % objective function
fun1and2 = @(x) [x(end), obj_fun(x)]; % all attributes (including parameter or parameter function)
par = [1];                          % the parameter function is in the first column of fun1and2
dom = [2];                          % the dominator function is in the second column of fun1and2
% Bound constraints on X = [hL, t, fd1, fd2, load]
lb = [0.05,  0.5, 0.1, 0.1,  2.5];
ub = [ 0.5,  2.5, 0.9, 0.9, 10.0];
nvars = length(lb);                 % the number of design variables
A = []; b = []; Aeq = []; beq = []; % linear inequality and equality constraints
Generations = 50;                   % maximum number of generations
PopulationSize = 50;                % maximum populations
ref_bnds = [2.5,10;300,3000];     % reference bounds for HVI measure
nonlcon = @(x) abaqus_constr(x,h0,L,yielding_stress,Disp_allowed); % nonlinear constraint function
classif_err_allowance = 0.2;
TolCon = 1e-6;

%% Section 4: Define 'options'
if W_constr_handling==2
    % Initialize population
    [X_new, feasible, initial, initialfeasible, cEval] = Initialize_LHS_GP(lb,ub,...
        A,b,Aeq,beq,nonlcon,PopulationSize,nvars,'DontPlot',classif_err_allowance,TolCon);
    % Define options
    options = p3gaoptimset('Generations',Generations,'PopulationSize',PopulationSize,...
        'ViewProgress',true,'InitialPopulation',X_new,...
        'CrossoverFraction',0.8,'MutationFraction',0.05);
    options.initial = initial;
    options.initialfeasible = initialfeasible;
    options.GenerationPlots = 1;
    options.classif_err_allowance = classif_err_allowance;
    options.cEval = cEval;
    options.TolCon = TolCon;
else
    % Define options
    options = p3gaoptimset('Generations',Generations,'PopulationSize',PopulationSize,...
        'ViewProgress',true,'CrossoverFraction',0.8,'MutationFraction',0.05);
end

% Save data by generation in a specified directionary
options.Log = true;
options.GenerationData = true;
if options.GenerationData
    if W_constr_handling==2
        options.SaveLoc = strcat('data\test1_w');
    elseif W_constr_handling==3
        options.SaveLoc = strcat('data\test3_w_CTAEA');
    elseif W_constr_handling==4
        options.SaveLoc = strcat('data\test2_w_cDPEA');
    else
        options.SaveLoc = strcat('data\test1_wo');
    end
    if ~exist(options.SaveLoc, 'dir')
        mkdir(options.SaveLoc)
    end
    options.GenerationDataIncrement = 1;
end
options.Dominance = 'hch'; % determine if you want to use HCH-based dominance
options.hvi_options.phi = 45; % Hypercone angle
options.ref_bnds = ref_bnds;

% if W_constr_handling
%     % Save the plot for initial population
%     figure(1);
%     title(['Initial Population'])
%     saveas(figure(1),strcat(options.SaveLoc,'\pop_sample',num2str(0),'.png'))
%     % close(gcf)
% end

%% Section 5: Run P3GA
% Minimize using P3GA
if W_constr_handling==3
    [x_opt,fval,M] = p3ga_TAEA(fun1and2,dom,par,nvars,A,b,Aeq,beq,lb,ub,nonlcon,options);
elseif W_constr_handling==4
    [x_opt,fval,M] = p3ga_cDPEA(fun1and2,dom,par,nvars,A,b,Aeq,beq,lb,ub,nonlcon,options);
else
    [x_opt,fval,M] = p3ga(fun1and2,dom,par,nvars,A,b,Aeq,beq,lb,ub,nonlcon,options);
end

figure(gcf); xlabel('Generation','Interpreter','Latex'); ylabel('pHVI','Interpreter','Latex');
saveas(figure(gcf),strcat(options.SaveLoc,'\p3ga_pHVI.fig'))
close(gcf)

figure(gcf); xlabel('$\theta$','Interpreter','Latex'); ylabel('$J$','Interpreter','Latex');
xlim([lb(end),ub(end)])
saveas(figure(gcf),strcat(options.SaveLoc,'\p3ga_result.fig'))
close(gcf)

pHVI_value = M.hvi(2,end);

%% Section 8: Save data
save(strcat(options.SaveLoc,'\matlab.mat'))

if W_constr_handling==2
    figure
    plot(1:Generations,M.exp_classif_err)
    xlabel('Generations'); ylabel('Expected value of classification errors')
end
