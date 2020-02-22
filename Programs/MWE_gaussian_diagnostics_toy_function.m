%
% Minimum working example to test Gaussian diagnostics idea
% using a custom build toy function
%
function MWE_gaussian_diagnostics_toy_function()

format short
close all

rng(202326) % initialize random number generator for reproducability

N = 2^(15);
f_c = randn(N, 1);
f_s = randn(N, 1);
f_0 = randn(1, 1);
  
% gaussian random function
function fval = f_rand(xpts,theta)

  kvec = (1:N)';
  argx = @(x) 2*pi*bsxfun(@times, kvec, x);
  f_c_ = @(x)(f_c./kvec).*cos(argx(x));
  f_s_ = @(x)(f_s./kvec).*sin(argx(x));
  f_ran = @(x,theta) prod((f_0 + theta * sum(f_c_(x) + f_s_(x) )), 2) ;
  [n,~]=size(xpts);
  fval=zeros(n,1);
  for i=1:n
    fval(i) = f_ran(xpts(i,:),theta);
  end
end

% bernPoly = @(x)(-x.*(1-x) + 1/6);

% fval = bernPoly(0.5);
% fval_ = bernoulli_series(0.5,2);

% fval = bernPoly(xlat_);
% fval_ = cubBayesLattice_g.bernoulli_series(xlat_,2);
% sum(sum(abs(fval - fval_)))


npts = 2^6;  % max 14
dim = 1;
shift = rand(1,dim);

[~,xlat_] = cubBayesLattice_g.simple_lattice_gen(npts,dim,shift,true);

fval = f_rand(xlat_, 1);
[H_,pValue_,W_] = swtest(fval)
if H_==1, Hval='false'; else, Hval='true'; end
fprintf('Shapiro-Wilk test: Normal=%s, pValue=%1.3e, W=%1.3f, n=%d\n', ...
  Hval, pValue_, W_, npts);

end
