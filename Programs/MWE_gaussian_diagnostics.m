%
% Minimum working example to test Gaussian diagnostics idea
%
function MWE_gaussian_diagnostics()

format short
close all

fNames = {'ExpCos','Keister','rand'};
ptransforms = {'C1','C1sin', 'none'};
fName = fNames{1};
ptransform = ptransforms{2};

npts = 2^4;  % max 14
dim = 1;
shift = rand(1,dim);

[~,xlat_] = simple_lattice_gen(npts,dim,shift,true);

if strcmp(fName,'ExpCos')
  integrand = @(x) exp(sum(cos(2*pi*x), 2));
elseif strcmp(fName, 'Keister')
  integrand = @(x) keisterFunc(x,dim,1/sqrt(2)); % a=0.8
else
  integrand = @(x) f_rand(x);
end

integrand_p = doPeriodTx(integrand, ptransform);

y = integrand_p(xlat_);
ftilde = fft(y);
ftilde(1) = 0;  % ftilde = \mV^H(\vf - m \vone)
if dim==1
  figure; scatter(xlat_, y, 10)
  title('Integrand')
end

thetaOpt = 1;  %0.9;
ftilde = real(ftilde);

  function create_plots(type)
    hFigNormplot = figure();
    set(hFigNormplot,'defaultaxesfontsize',16, ...
      'defaulttextfontsize',16, ... %make font larger
      'defaultLineLineWidth',0.75, 'defaultLineMarkerSize',8)
    if strcmp(type, 'normplot')
      normplot(w_ftilde)
    else
      qqplot(w_ftilde)
    end
    
    title(sprintf('%s n=%d Tx=%s r=%1.2f theta=%1.2f', ...
      fName, npts, ptransform, r, thetaOpt))
    saveas(hFigNormplot, sprintf('%s_%s_n-%d_Tx-%s_r-%d.png', ...
      type, fName, npts, ptransform, r))
  end

rVec = [2, ];

for r=rVec
  
  lambda = kernel(r, xlat_, thetaOpt);
  
  % apply transform
  % $\vZ = \frac 1n \mV \mLambda^{-\frac 12} \mV^H(\vf - m \vone)$
  % ifft also includes 1/n division
  w_ftilde = real(ifft(ftilde./sqrt(lambda)));
  
  % create_plots('normplot')
  create_plots('qqplot')
  
  % Shapiro-Wilk test
  [H, pValue, W] = swtest(w_ftilde);
  Hval='true';
  if H==true
    Hval='false';
  end
  fprintf('Shapiro-Wilk test: Normal=%s, pValue=%1.3f, W=%1.3f\n', Hval, pValue, W);
end

end


% gaussian random function
function fval = f_rand(xpts)
theta = 1;
rng(202326) % initialize random number generator for reproducability

N = 2^(15);
f_c = randn(N, 1);
f_s = randn(N, 1);
f_0 = randn(1, 1);
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

function Lambda = kernel(r, xun, theta)
constMult = 1;  % -(-1)^(r/2)*((2*pi)^r)/factorial(r);
if r == 2
  bernPoly = @(x)(-x.*(1-x) + 1/6);
elseif r == 4
  bernPoly = @(x)( ( (x.*(1-x)).^2 ) - 1/30);
else
  error('Bernoulli order=%d not implemented !', r);
end
kernelFunc = @(x) bernPoly(x);

temp_ = bsxfun(@times, (theta)*constMult, kernelFunc(xun));
C1 = prod(1 + temp_, 2);

% matlab's builtin fft is much faster and accurate
% eigenvalues must be real : Symmetric pos definite Kernel
Lambda = real(fft(C1));

end

function f = doPeriodTx(fInput, ptransform)

if strcmp(ptransform,'Baker')
  f=@(x) fInput(1-2*abs(x-1/2)); % Baker's transform
elseif strcmp(ptransform,'C0')
  f=@(x) fInput(3*x.^2-2*x.^3).*prod(6*x.*(1-x),2); % C^0 transform
elseif strcmp(ptransform,'C1')
  % C^1 transform
  f=@(x) fInput(x.^3.*(10-15*x+6*x.^2)).*prod(30*x.^2.*(1-x).^2,2);
elseif strcmp(ptransform,'C1sin')
  % Sidi C^1 transform
  f=@(x) fInput(x-sin(2*pi*x)/(2*pi)).*prod(2*sin(pi*x).^2,2);
elseif strcmp(ptransform,'C2sin')
  % Sidi C^2 transform
  psi3 = @(t) (8-9*cos(pi*t)+cos(3*pi*t))/16;
  psi3_1 = @(t) (9*sin(pi*t)*pi- sin(3*pi*t)*3*pi)/16;
  f=@(x) fInput(psi3(x)).*prod(psi3_1(x),2);
elseif strcmp(ptransform,'C3sin')
  % Sidi C^3 transform
  psi4 = @(t) (12*pi*t-8*sin(2*pi*t)+sin(4*pi*t))/(12*pi);
  psi4_1 = @(t) (12*pi-8*cos(2*pi*t)*2*pi+sin(4*pi*t)*4*pi)/(12*pi);
  f=@(x) fInput(psi4(x)).*prod(psi4_1(x),2);
elseif strcmp(ptransform,'none')
  % do nothing
  f=@(x) fInput(x);
else
  error('Error: Periodization transform %s not implemented', ptransform);
end

end


function [xlat,xpts_un,xlat_un,xpts] = simple_lattice_gen(n,d,shift,firstBatch)
gen_vec = [1, 433461, 315689, 441789, 501101, 146355, 88411, 215837, 273599 ...
  151719, 258185, 357967, 96407, 203741, 211709, 135719, 100779, ...
  85729, 14597, 94813, 422013, 484367]; %generator
z = gen_vec(1:d);

nmax = n;
nmin = 1 + n/2;
if firstBatch==true
  nmin = 1;
end
nelem=nmax-nmin+1;

if firstBatch==true
  brIndices=vdc(nelem)';
  xpts_un=mod(bsxfun(@times,(0:1/n:1-1/n)',z),1); % unshifted in direct order
else
  brIndices=vdc(nelem)'+1/(2*(nmin-1));
  xpts_un=mod(bsxfun(@times,(1/n:2/n:1-1/n)',z),1); % unshifted in direct order
  
end
xpts = mod(bsxfun(@plus,xpts_un,shift),1);  % shifted in direct order

xlat_un = mod(bsxfun(@times,brIndices',z),1);  % unshifted
xlat = mod(bsxfun(@plus,xlat_un,shift),1);  % shifted in VDC order
end

% van der Corput sequence in base 2
function q = vdc(n)
if n>1
  k=log2(n); % We compute the VDC seq part by part of power 2 size
  q=zeros(2^k,1);
  for l=0:k-1
    nl=2^l;
    kk=2^(k-l-1);
    ptind=repmat([false(nl,1);true(nl,1)],kk,1);
    q(ptind)=q(ptind)+1/2^(l+1);
  end
else
  q=0;
end
end

