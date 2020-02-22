%
% Minimum working example to test Gaussian diagnostics idea
%
function MWE_gaussian_diagnostics()

format short
close all

fNames = {'ExpCos','Keister','rand'};
ptransforms = {'C1','C1sin', 'none'};
fName = fNames{1};
ptransform = ptransforms{3}; 
whKer = 'Bern';  %'Exp'; %

npts = 2^6;  % max 14
dim = 1;
shift = rand(1,dim);

[~,xlat_] = cubBayesLattice_g.simple_lattice_gen(npts,dim,shift,true);

if strcmp(fName,'ExpCos')
  integrand = @(x) exp(sum(cos(2*pi*x), 2));
elseif strcmp(fName, 'Keister')
  integrand = @(x) keisterFunc(x,dim,1/sqrt(2)); % a=0.8
else
  integrand = @(x) f_rand(x, 1);
end
      
% integrand = @(x) prod(bernPoly(x), 2);

inputArgs = {'f',integrand,'fName',fName,'dim',dim, 'order',1, ....
  'ptransform',ptransform, 'stopAtTol',true, 'stopCriterion','MLE', ...
  'arbMean',true, 'visiblePlot',true};

obj=cubBayesLattice_g(inputArgs{:});

integrand_p = cubBayesLattice_g.doPeriodTx(integrand, ptransform);

y = integrand_p(xlat_);
ftilde = fft(y);
ftilde(1) = 0;  % ftilde = \mV^H(\vf - m \vone)
if dim==1
  figure; scatter(xlat_, y, 10)
  title('Integrand')
end

avoidCancelError = false;
thetaOpt = 0.9;
ftilde = real(ftilde);

rVec = [obj.order];
for r=rVec
  debugEnable = false;

  [lambda] = cubBayesLattice_g.dKernel(xlat_,obj.order,thetaOpt,avoidCancelError, ...
           obj.kernType,false,debugEnable);
  
  % apply transform
  % $\vZ = \frac 1n \mV \mLambda^{-\frac 12} \mV^H(\vf - m \vone)$
  w_ftilde = real(ifft(ftilde./sqrt(abs(real(lambda)))));

  hFigNormplot = figure();
  set(hFigNormplot,'defaultaxesfontsize',16, ...
    'defaulttextfontsize',16, ... %make font larger
    'defaultLineLineWidth',0.75, 'defaultLineMarkerSize',8)
  normplot(w_ftilde)
  
  plot_title = sprintf('%s n=%d Tx=%s r=%1.2f theta=%1.2f', ...
    fName, npts, ptransform, r, thetaOpt);
  title(plot_title)
  saveas(hFigNormplot, sprintf('normplot_%s_n-%d_Tx-%s_r-%d.png', ...
    fName, npts, ptransform, r))
  
  % Shapiro-Wilk test
  [H, pValue, W] = swtest(w_ftilde);
  Hval='true';
  if H==true
    Hval='false';
  end
  fprintf('Shapiro-Wilk test: Normal=%s, pValue=%1.3f, W=%1.3f\n', Hval, pValue, W);
end

end




