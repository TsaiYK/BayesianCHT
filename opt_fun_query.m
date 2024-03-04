function [xQ,fval,eval_func] = opt_fun_query(fun,x0,lb,ub)
    % This function is to find the query for classification by GA optimizer
    % The objective function is to maximize the entropy by BALD, Houlsby,
    % N., Huszár, F., Ghahramani, Z., & Lengyel, M. (2011). Bayesian active 
    % learning for classification and preference learning. arXiv preprint arXiv:1112.5745.
    % 
    % Input:
    % fun: objective function for the optimization = - (predictive output variance
    % for classification) because the optimal query would be the point that
    % have the highest uncertainty for the classification
    % x0: initial for optimization
    A = []; b = []; Aeq = []; beq = [];
    nonlcon = [];
    options = optimoptions('ga','FunctionTolerance',1e-3,'MaxGenerations',min(max(10*(length(x0)-1),10),50),...
        'PopulationSize',min(max(10*(length(x0)-1),10),50),'Display','off');
    [xQ,fval,exitflag,output] = ga(fun,length(lb),A,b,Aeq,beq,lb,ub,nonlcon,options);
    eval_func = output.funccount;
    
    if exitflag~=1 && exitflag~=0
        fval = [];
        xQ = [];
    end
end