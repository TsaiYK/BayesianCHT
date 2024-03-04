function [ymu, ys2, fmu, fs2, lp] = gp_Eval(hyp, inf, mean, cov, lik, post, x, y, xs, ys)

%%
if isempty(mean), mean = {@meanZero}; end                     % set default mean
if ischar(mean) || isa(mean, 'function_handle'), mean = {mean}; end  % make cell
if isempty(cov), error('Covariance function cannot be empty'); end  % no default
if ischar(cov) || isa(cov,'function_handle'), cov  = {cov};  end     % make cell
cstr = cov{1}; if isa(cstr,'function_handle'), cstr = func2str(cstr); end
if (strcmp(cstr,'covFITC') || strcmp(cstr,'apxSparse')) && isfield(hyp,'xu')
  cov{3} = hyp.xu;                                                   %use hyp.xu
end
if isempty(inf), inf = {@infGaussLik}; end        % set default inference method
if ischar(inf),  inf = str2func(inf);  end        % convert into function handle
if ischar(inf) || isa(inf,'function_handle'), inf = {inf};  end      % make cell
istr = inf{1}; if isa(istr,'function_handle'), istr = func2str(istr); end
if strcmp(istr,'infPrior')
  istr = inf{2}; if isa(istr,'function_handle'), istr = func2str(istr); end
end
if isempty(lik), lik = {@likGauss}; end                        % set default lik
if ischar(lik) || isa(lik,'function_handle'), lik = {lik};  end      % make cell
lstr = lik{1}; if isa(lstr,'function_handle'), lstr = func2str(lstr); end

%%
if ~isfield(hyp,'mean'), hyp.mean = []; end        % check the hyp specification
% if eval(feval(mean{:})) ~= numel(hyp.mean)
%   error('Number of mean function hyperparameters disagree with mean function')
% end
if ~isfield(hyp,'cov'), hyp.cov = []; end
% if eval(feval(cov{:})) ~= numel(hyp.cov)
%   error('Number of cov function hyperparameters disagree with cov function')
% end
if ~isfield(hyp,'lik'), hyp.lik = []; end
% if eval(feval(lik{:})) ~= numel(hyp.lik)
%   error('Number of lik function hyperparameters disagree with lik function')
% end

%%
alpha = post.alpha; L = post.L; sW = post.sW;
  if issparse(alpha)                  % handle things for sparse representations
    nz = alpha ~= 0;                                 % determine nonzero indices
    if issparse(L), L = full(L(nz,nz)); end      % convert L and sW if necessary
    if issparse(sW), sW = full(sW(nz)); end
  else nz = true(size(alpha,1),1); end               % non-sparse representation
  if isempty(L)                       % in case L is not provided, we compute it
    K = feval(cov{:}, hyp.cov, x(nz,:));
    L = chol(eye(sum(nz))+sW*sW'.*K);
  end
  %verify whether L contains valid Cholesky decomposition or something different
  Lchol = isnumeric(L) && all(all(tril(L,-1)==0)&diag(L)'>0&isreal(diag(L))');
  ns = size(xs,1);                                       % number of data points
  if strncmp(cstr,'apxGrid',7), xs = apxGrid('idx2dat',cov{3},xs); end  % expand
  nperbatch = 1000;                       % number of data points per mini batch
  nact = 0;                       % number of already processed test data points
  ymu = zeros(ns,1); ys2 = ymu; fmu = ymu; fs2 = ymu; lp = ymu;   % allocate mem
  while nact<ns               % process minibatches of test cases to save memory
    id = (nact+1):min(nact+nperbatch,ns);               % data points to process
    kss = feval(cov{:}, hyp.cov, xs(id,:), 'diag');              % self-variance
    if strcmp(cstr,'covFITC') || strcmp(cstr,'apxSparse')    % cross-covariances
      Ks = feval(cov{:}, hyp.cov, x, xs(id,:)); Ks = Ks(nz,:); % res indep. of x
    else
      Ks = feval(cov{:}, hyp.cov, x(nz,:), xs(id,:));        % avoid computation
    end
    ms = feval(mean{:}, hyp.mean, xs(id,:));
    N = size(alpha,2);  % number of alphas (usually 1; more in case of sampling)
    Fmu = repmat(ms,1,N) + Ks'*full(alpha(nz,:));        % conditional mean fs|f
    fmu(id) = sum(Fmu,2)/N;                                   % predictive means
    if Lchol    % L contains chol decomp => use Cholesky parameters (alpha,sW,L)
      V  = L'\(repmat(sW,1,length(id)).*Ks);
      fs2(id) = kss - sum(V.*V,1)';                       % predictive variances
    else                % L is not triangular => use alternative parametrisation
      if isnumeric(L), LKs = L*Ks; else LKs = L(Ks); end    % matrix or callback
      fs2(id) = kss + sum(Ks.*LKs,1)';                    % predictive variances
    end
    fs2(id) = max(fs2(id),0);   % remove numerical noise i.e. negative variances
    Fs2 = repmat(fs2(id),1,N);     % we have multiple values in case of sampling
    if nargin<9
      [Lp, Ymu, Ys2] = feval(lik{:},hyp.lik,[],   Fmu(:),Fs2(:));
    else
      Ys = repmat(ys(id),1,N);
      [Lp, Ymu, Ys2] = feval(lik{:},hyp.lik,Ys(:),Fmu(:),Fs2(:));
    end
    lp(id)  = sum(reshape(Lp, [],N),2)/N;    % log probability; sample averaging
    ymu(id) = sum(reshape(Ymu,[],N),2)/N;          % predictive mean ys|y and ..
    ys2(id) = sum(reshape(Ys2,[],N),2)/N;                          % .. variance
    nact = id(end);          % set counter to index of last processed data point
  end
  if nargin<9
    varargout = {ymu, ys2, fmu, fs2, [], post};        % assign output arguments
  else
    varargout = {ymu, ys2, fmu, fs2, lp, post};
  end