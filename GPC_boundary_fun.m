function f = GPC_boundary_fun(x, X, y, hyp, inffunc, meanfunc, covfunc, likfunc, post)
    [ymu, ys2, fmu, fs2, lp] = gp_Eval(hyp, inffunc, meanfunc, covfunc, likfunc, post, X, y, x, 1);
    f = exp(lp);
end