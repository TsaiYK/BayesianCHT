function [X_new, feasible_repaired, X, feasible, cEval] = Initialize_LHS_GP(lb,ub,A,b,Aeq,beq,nonlcon,Npop,nvars,PlotOrNot,classif_err_allowance,TolCon,X_initial)
% This function is to generate the initial populations using LHS and GP 
% techniques. The goal is to generate the initial population for P3GA 
% through identifying the nonlinear constraint boundary and repairing those
% infeasible points to the predicted boundary. 

delta_x = ub-lb;
GPC_success = false;
counter = 1; X_train = []; y_train = [];
while ~GPC_success
    if isempty(X_initial)
        % Latin Hypercube Sampling
        X{counter} = lhsdesign(Npop,nvars).*repmat(delta_x,Npop,1)+repmat(lb,Npop,1);
    else
        X{counter} = X_initial;
    end

    % initialize the feasible array
    feasible = true([Npop,1]);

    % evaluate the nonlinear constraint function
    cEval = 0;
    for i = 1:Npop
        feasible(i,:) = isfeasible(X{counter}(i,:),A,b,Aeq,beq,lb,ub,nonlcon,TolCon);
        cEval = cEval+1;
    end
    
    %% Classify based on the feasibility
    cc = feasible==false; % the ones violate the constraints (c>0)
    infeasible_index = find(sum(cc,2)>0);
    if size(infeasible_index,2) == 1
        infeasible_index = infeasible_index'; % 'row' should be a row array
    end

    [X_feasible,X_infeasible] = feasible_infeasible_plot(X{counter},feasible);
    if isempty(X_infeasible) || length(X_infeasible)/Npop<0.1
        X_new = X_feasible;
        feasible_repaired = feasible;
        break;
    else
        %% GP + shift
        y{counter} = ones(size(feasible));
        y{counter}(feasible==false) = -1; % the gp model only accepts 1 or -1
        for i = 1:counter
            X_train = [X_train;X{counter}];
            y_train = [y_train;y{counter}];
        end
        [GP_fun, Boundary_fun, hyp, inffunc, meanfunc, covfunc, likfunc, post] = train_GPC(X_train,y_train);

        % Calculate the expected classification error
        for j = 1:size(X{counter},1)
            [~, ~, ~, ~, lp] = gp_Eval(hyp, inffunc, meanfunc, covfunc, likfunc, post, X_train, y_train, X{counter}(j,:), 1);
            prob(j) = exp(lp);
            if sum(j==infeasible_index)==1 % should be 0
                classif_err(j) = prob(j);
            else % should be 1
                classif_err(j) = 1-prob(j);
            end
        end
        exp_classif_err = mean(classif_err);

        if isempty(X_infeasible)
            NoNeedShift = true;
        else
            NoNeedShift = false;
        end

        if ~NoNeedShift
            [X_infeasible_shift,~,shift_index] = RepairInfeasible...
                (GP_fun,Boundary_fun,X_infeasible,X_train,y_train,lb,ub,...
                hyp, inffunc, meanfunc, covfunc, likfunc, post, exp_classif_err, classif_err_allowance);
            if isempty(X_infeasible_shift)
                GPC_success = false;
                counter = counter+1;
                continue
            else
                GPC_success = true;
            end

            feasible_repaired = true([1,Npop]);
            feasible_repaired = feasible; X_new = X{counter};
            k1 = 1; k2 = 1;
            for i = infeasible_index
                if sum(k1==shift_index)==1 % this point was shifted
                    feasible_repaired(i,:) = isfeasible(X_infeasible_shift(k2,:),A,b,Aeq,beq,lb,ub,nonlcon,TolCon);
                    X_new(i,:) = X_infeasible_shift(k2,:);
                    cEval = cEval+1;
                    k2 = k2+1;
                end
                k1 = k1+1;
            end
            cc_shifted = feasible_repaired==false; % the ones violate the constraints (c>0)
            infeasible_shift_index = sum(cc_shifted,2)>0;
            feasible_repaired(infeasible_shift_index) = false;
        end
        counter = counter+1;
        if counter>5
            break
        end
    end
end
if counter>5
    warning('No feasible found for the initial populations!')
end

%% Plots
if strcmp(PlotOrNot,'Plot')==1 && length(lb)==2
    plotGP(X_feasible,X_infeasible,lb,ub,X_train,y_train, hyp, inffunc, meanfunc, covfunc, likfunc, post); % visualize
end