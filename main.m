% This script is to define the problem, CHT, and any parameters for running 
% the algorithm and saving the results and call the functions to run P3GA. 
clear
clc
close all

set(0,'DefaultTextInterpreter','latex'); % change the text interpreter
set(0,'DefaultLegendInterpreter','latex'); % change the legend interpreter
set(0,'DefaultAxesTickLabelInterpreter','latex'); % change the tick interpreter
set(0, 'DefaultLineLineWidth', 1);
set(groot,'defaultAxesXGrid','on')
set(groot,'defaultAxesYGrid','on')
set(0,'defaultAxesFontSize',14)

addpath(genpath(pwd));

%% Section 1: Setting (define functions given problem type and others)
% Define scalability:
% m: number of objective functions
% p: number of parameter functions

classif_err_allowance = 0.2;
m = 2; p = 0;
n = m+p;
ProbType = 3;
W_constr_handling = 3;
run_trial = 1;

%% Section 2: Define Objective, Parameter, and Constraint Functions
obj_fun = @(x) objectivefunA(x,m,p);            % objective function
para_fun = @(x) parameterfunA(x,p);             % parameter function
Prob.fun1and2 = @(x) [para_fun(x), obj_fun(x)]; % all attributes (including parameter or parameter function)
Prob.par = [1:p];                               % index for the parameter functions in fun1and2
Prob.dom = [(p+1):n];                           % index for the dominance functions in fun1and2
Prob.lb = zeros(1,n);                           % lower bound
Prob.ub = ones(1,n);                            % upper bound
Prob.nvars = length(Prob.lb);                   % the number of design variables
Prob.Generations = 30;                          % maximum number of generations
Prob.PopulationSize = 100;                      % maximum populations

% Determine nonlinear constraint function given problem type
switch ProbType
    case 1
        Prob.nonlcon = @(x) nonlconfunA1(x,m,p);
    case 2
        Prob.nonlcon = @(x) nonlconfunA2(x,m,p);
    case 3
        Prob.nonlcon = @(x) nonlconfunA3(x,m,p);
    case 4
        Prob.nonlcon = @(x) nonlconfunA4(x,m,p);
end

% Reference point = [x_lb,x_ub;y_lb,y_ub]
if n==2
    Prob.ref_bnds = [0,1;6,16];
elseif n==3
    Prob.ref_bnds = [0,1;0,1;9,25];
elseif n==4
    Prob.ref_bnds = [0,1;0,1;0,1;13,33];
end

%% Section 3: Add paths
if W_constr_handling==3
    rmpath(genpath('PlatEMO\PlatEMO\Algorithms\Multi-objective optimization\c-DPEA'))
    addpath(genpath('PlatEMO\PlatEMO\Algorithms\Multi-objective optimization\C-TAEA'))
    addpath(genpath('PlatEMO\PlatEMO\Algorithms\Utility functions'))
    addpath(genpath('PlatEMO\PlatEMO\Problems'))
elseif W_constr_handling==4
    rmpath(genpath('PlatEMO\PlatEMO\Algorithms\Multi-objective optimization\C-TAEA'))
    addpath(genpath('PlatEMO\PlatEMO\Algorithms\Multi-objective optimization\c-DPEA'))
    addpath(genpath('PlatEMO\PlatEMO\Algorithms\Utility functions'))
    addpath(genpath('PlatEMO\PlatEMO\Problems'))
end
 
%% Section 4: Run
[pHVI_value,x_opt,fval,M,exp_classif_err] = run_p3ga(Prob,W_constr_handling,classif_err_allowance,ProbType,run_trial);



