function [xP,fval,eval_func] = opt_fun_project(x0,lb,ub,funP)
    % This function is to find the nearest point on the predicted boundary
    % Input:
    % funP: objective function for minimizing the distance between the
    % infeasible point and the predicted boundary
    % x0: initial for optimization
    delta = ub-lb;
    xm = lb+delta/2;
    x0 = rand(1,length(ub)).*delta+(xm-delta/2);
    fun = @(x) obj_fun(x,x0);
    A = []; b = []; Aeq = []; beq = [];
    nonlcon = @(x) nonlcon_GPC(x, funP);
    options = optimoptions('fmincon','Display','off','MaxFunctionEvaluations',...
        100*(length(x0)-1),'ConstraintTolerance',1e-3);
    [xP,fval,exitflag,output] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
    eval_func = output.funcCount;
    
    if exitflag<0
        fval = [];
        xP = [];
    end
    
    function f = obj_fun(x,x0)
        f = norm(x-x0);
    end
    function [c,ceq] = nonlcon_GPC(x, funP)
        prob = funP(x);
        c = 0.5-prob;
        ceq = [];
    end
end