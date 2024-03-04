function [xQ,xQ_infea,shift_index] = RepairInfeasible(funQ,funP,x2,X,y,lb,ub,...
    hyp, inffunc, meanfunc, covfunc, likfunc, post, exp_classif_err, classif_err_allowance)
% The function is to reapir the infeasible points using Bayesian methods 
% (Bayesian active learning and Gaussian process classifier). The repair
% method depends on if exp_classif_err<classif_err_allowance. If so, use 
% `opt_fun_project' (move the infeasible points to the nearest constraint 
% boundaries). If not, use `opt_fun_query' to query new points that are
% the most informative and optimally improve the model accuracy
% Input: 
% funQ: the entropy function entropy_fun, defined by BALD . This function uses the trained GPC model.
% funP: the lost function, GPC_boundary_fun, to be minimized. This function is to find the closest point on the constraint boundary.
% x2: infeasible points to be repaired
% Output:
% xQ: repaired points
% xQ_infea: unrepaired points (for some reasons, the points fail to be repaired)

% Count the number of infeasible solutions
num_infeasible = size(x2,1);
fprintf('Number of infeasible points: %d\n',num_infeasible);

% Initialize X_new
X_new = X; 

% Initialize the infeasible points that fail repairing
xQ_infea = []; shift_index = [];

% Determine quering or projecting based on the expected classif error
if exp_classif_err<classif_err_allowance
    disp('Projecting...')
    num = 0;
    for i = 1:num_infeasible
        % projecting to the nearest feasible point
        [xQ_tmp,~,~] = opt_fun_project(x2(i,:),lb,ub,funP);
        
        if isempty(xQ_tmp) % unsuccessful projecting
            xQ_infea = [xQ_infea;x2(i,:)];
        else
            X_new = [X_new;xQ_tmp];
            shift_index = [shift_index,i];
            num = num+1;
        end
    end
    fprintf('Number of moved points: %d\n',num);
else % keep querying
    disp('Querying...')
    num = 0;
    for i = 1:num_infeasible
        % Optimization (to find the point that has the highest entropy)
        [xQ_tmp,~,~,~] = opt_fun_query(funQ,x2(i,:),lb,ub);
        
        % Evaluate GPC model
        [~, ~, ~, ~, lp] = gp_Eval(hyp, inffunc, meanfunc, covfunc, likfunc, post, X, y, xQ_tmp, 1);
        
        % probability of y=+1
        p = exp(lp);
        
%         if abs(p-0.5)>0.1
%             xQ_infea = [xQ_infea;x2(i,:)];
        if i==1 && isempty(xQ_tmp) % the first infeasible and unsuccessful query
            xQ_infea = [xQ_infea;x2(i,:)];
        else % not the first infeasible, we need to determine the closeness
            X_new = [X_new;xQ_tmp];
            shift_index = [shift_index,i];
            num = num+1;
        end
    end
    fprintf('Number of queried points: %d\n',num);
end
n_train = size(X,1);
xQ = X_new(n_train+1:end,:);
% fval = min(fval);
