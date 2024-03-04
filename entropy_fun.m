function f = entropy_fun(x, X, y, hyp, inffunc, meanfunc, covfunc, likfunc, post)
% Funtion of entropy used as the objective function for BALD active learning
% Reference:
% N. Houlsby, F. Husz´ar, Z. Ghahramani, and M. Lengyel, “Bayesian active
% learning for classification and preference learning,” arXiv preprint
% arXiv:1112.5745, 2011.
    [ymu, ys2, fmu, fs2, lp] = gp_Eval(hyp, inffunc, meanfunc, covfunc, likfunc, post, X, y, x, 1);
    p = exp(lp);
    h = -p*log10(p)-(1-p)*log10(1-p);
    C = sqrt(pi*log(2)/2);
    E_fH = C/sqrt(ys2+C^2)*exp(-ymu^2/2/(ys2+C^2));
    f = -(h-E_fH); % to maximize the entropy
end