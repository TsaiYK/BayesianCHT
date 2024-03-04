function [fun1, fun2, hyp, inffunc, meanfunc, covfunc, likfunc, post] = train_GPC(X,y)
% This function is to train the gaussian process classification (GPC) model 
% using the toolkit developed by Rasmussen and Williams . In the function,
% the inference function, mean function, etc. are defined and use the dataset
% (X,y) to produce a GP classifier. Then, the loss functions for repair 
% mechanism are defined (maximize entropy function and minimize the distance 
% to the approximate boundary).
% 
% Input: 
% X: (N-by-n) array, where N is the number of samples and n is the dimension of the design variables
% y: (N-by-1) array with +1 or -1
% Output:
% fun1
% fun2
% GP model outputs (hyp, inffunc, meanfunc, covfunc, likfunc, post)

    inffunc = @infLaplace;
    meanfunc = @meanConst; hyp.mean = 0;
    covfunc = @covSEard; ell = 1.0; sf = 1.0; hyp.cov = log([ell*ones(1,size(X,2)) sf]);
    likfunc = @likErf;
    
    hyp = minimize(hyp, @gp, -40, inffunc, meanfunc, covfunc, likfunc, X, y);
    % [ymu, ys2, fmu, fs2, lp] = gp(hyp, inffunc, meanfunc, covfunc, likfunc, X, y, t, ones(n, 1));
    % [ymu, ys2, fmu, fs2, lp] = gp(hyp, inffunc, meanfunc, covfunc, likfunc, X, y, [-2.4,-0.6], ones(1, 1));
    [nlZ, dnlZ, post] = gp(hyp, inffunc, meanfunc, covfunc, likfunc, X, y);
    %   ymu      column vector (of length ns) of predictive output means
    %   ys2      column vector (of length ns) of predictive output variances
    %   fmu      column vector (of length ns) of predictive latent means
    %   fs2      column vector (of length ns) of predictive latent variances
    %   lp       column vector (of length ns) of log predictive probabilities
    
    fun1 = @(x) entropy_fun(x, X, y, hyp, inffunc, meanfunc, covfunc, likfunc, post); % min -entroy
    fun2 = @(x) GPC_boundary_fun(x, X, y, hyp, inffunc, meanfunc, covfunc, likfunc, post); % find boundary
end