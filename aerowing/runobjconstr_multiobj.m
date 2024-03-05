function [x_opt,fval,M] = runobjconstr_multiobj(x0,lb,ub,W_constr_handling)
xLast = []; % Last place computeall was called
myf1 = []; % Use for objective at xLast
myf2 = []; % Use for objective at xLast
myc = []; % Use for nonlinear inequality constraint
myceq = []; % Use for nonlinear equality constraint

% Data = load('C:\Users\yktsai0121\Desktop\AERO604_abaqus\CaseStudy_AeroWing\data\p3ga_results_w_CHT\initial_population.mat');
% X_new = Data.X_new;
% feasible = Data.feasible;
% initial = Data.initial;
% initialfeasible = Data.initialfeasible;
% cEval = Data.cEval;

A = []; b = []; Aeq = []; beq = [];
TolCon = 1e-6; classif_err_allowance = 0.2;
nvars = length(lb);

Generations = 20;
PopulationSize = 30;
% SaveLoc_dir = 'data\p3ga_results_w_CHT6';
if W_constr_handling==2
    SaveLoc_dir = strcat('data\p3ga_results_w_CHT6');
elseif W_constr_handling==3
    SaveLoc_dir = strcat('data\p3ga_results_w_CTAEA1');
elseif W_constr_handling==4
    SaveLoc_dir = strcat('data\p3ga_results_w_cDPEA1');
else
    SaveLoc_dir = strcat('data\p3ga_results_wo_CHT6');
end
    
fun = @(x) [objfun1(x),objfun2(x),x(8)]; % the attribute function: obj1, obj2, para
nonlcon = @(x) constr(x); % the constraint function, nested below

if W_constr_handling~=2
    options = p3gaoptimset('Generations',Generations,'PopulationSize',PopulationSize,...
        'ViewProgress',true);
    options.SaveLoc = SaveLoc_dir;
else
    % Initialize population
    [X_new, feasible, initial, initialfeasible, cEval] = Initialize_LHS_GP(lb,ub,...
            A,b,Aeq,beq,nonlcon,PopulationSize,nvars,'DontPlot',classif_err_allowance,TolCon,[]);
    save(strcat(SaveLoc_dir,'\initial_population.mat'))
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
    options.SaveLoc = SaveLoc_dir;
end
options.Log = true;
options.GenerationData = true;
options.GenerationDataIncrement = 1;
options.Dominance = 'hch'; % determine if you want to use HCH-based dominance
options.hvi_options.phi = 90; % Hypercone angle
options.ref_bnds = [20,150;-60,0;40,65];

% Call fmincon
% [x_opt,fval,Result.exitflag,Result.output] = fmincon(fun,x0,[],[],[],[],lb,ub,cfun,opts);
% [x_opt,fval,Result.exitflag,Result.output,Result.population,Result.scores] = ga(fun,length(x0),[],[],[],[],lb,ub,cfun,opts);

[x_opt,fval,M] = p3ga(fun,[1,2],[3],length(x0),[],[],[],[],lb,ub,nonlcon,options);
figure(1)
xlabel('Objective 1: Mass [kg]')
ylabel('Objective 2: -Safety factor')
zlabel('Parameter: air speed [m/s]')
figure(1)
saveas(figure(gcf),strcat(options.SaveLoc,'\p3ga_result.fig'))
saveas(figure(gcf),strcat(options.SaveLoc,'\p3ga_result.png'))



    function y = objfun1(x)
        xDesign = x;
        xDesign(1) = round(x(1));
        xDesign(7) = round(x(7));
        if ~isequal(xDesign,xLast) % Check if computation is necessary
%             [myf,myc,myceq] = computeall_weighted_obj(xDesign);
            [myf1,myf2,myc,myceq,OutputVal] = computeall_multiobj(xDesign);
            xLast = xDesign;
        end
        % Now compute objective function
        y = myf1;
    end

    function y = objfun2(x)
        xDesign = x;
        xDesign(1) = round(x(1));
        xDesign(7) = round(x(7));
        if ~isequal(xDesign,xLast) % Check if computation is necessary
%             [myf,myc,myceq] = computeall_weighted_obj(xDesign);
            [myf1,myf2,myc,myceq,OutputVal] = computeall_multiobj(xDesign);
            xLast = xDesign;
        end
        % Now compute objective function
        y = myf2;
    end

    function [c,ceq] = constr(x)
        xDesign = x;
        xDesign(1) = round(x(1));
        xDesign(7) = round(x(7));
        if ~isequal(xDesign,xLast) % Check if computation is necessary
%             [myf,myc,myceq] = computeall_weighted_obj(xDesign);
            [myf1,myf2,myc,myceq,OutputVal] = computeall_multiobj(xDesign);
            xLast = xDesign;
        end
        % Now compute constraint functions
        c = myc; % In this case, the computation is trivial
        ceq = myceq;
    end

end