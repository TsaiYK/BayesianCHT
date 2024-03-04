function [pHVI_value,x_opt,fval,M,exp_classif_err] = run_p3ga(Prob,W_constr_handling,classif_err_allowance,ProbType,run_trial)
% This function is to run P3GA given the defined input arguments
% (Prob,W_constr_handling,classif_err_allowance,ProbType,run_trial)

lb = Prob.lb; ub = Prob.ub; fun1and2 = Prob.fun1and2;
dom = Prob.dom; par = Prob.par; nvars = Prob.nvars; nonlcon = Prob.nonlcon;
Generations = Prob.Generations; PopulationSize = Prob.PopulationSize;
ref_bnds = Prob.ref_bnds; p = length(par);
TolCon = 1e-6;
% Linear constraints
A = []; b = []; Aeq = []; beq = [];

%% Section 1: Define 'options'
if W_constr_handling==2
    % Initialize population
    [X_new, feasible, initial, initialfeasible, cEval] = Initialize_LHS_GP(lb,ub,...
        A,b,Aeq,beq,nonlcon,PopulationSize,nvars,'DontPlot',classif_err_allowance,TolCon,[]);
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
        options.SaveLoc = strcat('data\A\w_constr_multipleA_Type',...
            num2str(ProbType),'_n',num2str(nvars),'p',num2str(p),'\trial',num2str(run_trial));
    elseif W_constr_handling==3
        options.SaveLoc = strcat('data\A\wTAEA_constr_multipleB_Type',...
            num2str(ProbType),'_n',num2str(nvars),'p',num2str(p),'\trial',num2str(run_trial));
    elseif W_constr_handling==4
        options.SaveLoc = strcat('data\A\wDPEA_constr_multipleA_Type',...
            num2str(ProbType),'_n',num2str(nvars),'p',num2str(p),'\trial',num2str(run_trial));
    else
        options.SaveLoc = strcat('data\A\wo_constr_multipleA_Type',...
            num2str(ProbType),'_n',num2str(nvars),'p',num2str(p),'\trial',num2str(run_trial));
    end
    if ~exist(options.SaveLoc, 'dir') % if the folder does not exist, create it
        mkdir(options.SaveLoc)
    end
    options.GenerationDataIncrement = 5;
end
options.Dominance = 'hch'; % determine if you want to use HCH-based dominance
options.hvi_options.phi = 45; % Hypercone angle
options.ref_bnds = ref_bnds;


%% Section 2: Run P3GA
% Minimize using P3GA
if W_constr_handling==3
    [x_opt,fval,M] = p3ga_TAEA(fun1and2,dom,par,nvars,A,b,Aeq,beq,lb,ub,nonlcon,options);
elseif W_constr_handling==4
    [x_opt,fval,M] = p3ga_cDPEA(fun1and2,dom,par,nvars,A,b,Aeq,beq,lb,ub,nonlcon,options);
elseif W_constr_handling==2
    [x_opt,fval,M] = p3ga_Bayesian(fun1and2,dom,par,nvars,A,b,Aeq,beq,lb,ub,nonlcon,options);
else
    [x_opt,fval,M] = p3ga(fun1and2,dom,par,nvars,A,b,Aeq,beq,lb,ub,nonlcon,options);
end
figure(gcf); xlabel('Generation','Interpreter','Latex'); ylabel('pHVI','Interpreter','Latex');
saveas(figure(gcf),strcat(options.SaveLoc,'\p3ga_pHVI.fig'))
close(gcf)

if nvars<=3 % P3GA generates nondominated solutions if nvars<=3
    figure(gcf); xlabel('$\theta$','Interpreter','Latex'); ylabel('$J$','Interpreter','Latex');
    xlim([lb(1),ub(1)])
    saveas(figure(gcf),strcat(options.SaveLoc,'\p3ga_result.fig'))
    close(gcf)
end

pHVI_value = M.hvi(2,end);

%% Section 3: Save data
if W_constr_handling==2
    exp_classif_err = M.exp_classif_err;
    figure
    plot(1:Generations,M.exp_classif_err)
    xlabel('Generations'); ylabel('Expected value of classification errors')
    saveas(figure(gcf),strcat(options.SaveLoc,'\p3ga_expected_classif_err.fig'))
    close(gcf)
else
    exp_classif_err = [];
end

save(strcat(options.SaveLoc,'\matlab.mat'))
