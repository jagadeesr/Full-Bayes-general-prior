%
% Minimum working example to test Gaussian diagnostics idea
%
function MWE_gaussian_diagnostics()

format short
close all
gail.InitializeDisplay

whEx = 3;
fNames = {'ExpCos','Keister','rand'};
ptransforms = {'C1','C1sin', 'none'};
fName = fNames{whEx};
ptransform = ptransforms{whEx};

npts = 2^4;  % max 14
dim = 3;
rVec = 1.75;
%rVec = [1.75 2.45 3.15]; %vector of possible r values

for r=rVec
  
  %parameters for random function
  %seed = 202326;
  seed = randi([1,1e6],1,1)
  rfun = r/2;
  f_mean = 0;
  f_std_a = 1;
  f_std_b = 1;
  theta = (f_std_a/f_std_b)^2;
  
  shift = rand(1,dim);
  
  [~,xlat] = simple_lattice_gen(npts,dim,shift,true);
  
  if strcmp(fName,'ExpCos')
    integrand = @(x) exp(sum(cos(2*pi*x), 2));
  elseif strcmp(fName, 'Keister')
    integrand = @(x) keisterFunc(x,dim,1/sqrt(2)); % a=0.8
  elseif strcmp(fName, 'rand')
    integrand = @(x) f_rand(x, rfun, f_std_a, f_std_b, f_mean, seed);
  else
    error('Invalid function name')
  end
  
  integrand_p = doPeriodTx(integrand, ptransform);
  
  y = integrand_p(xlat); %function data
  ftilde = fft(y); %fourier coefficients
  ftilde(1) = 0;  % ftilde = \mV^H(\vf - m \vone), subtract mean
  if dim==1
    hFigIntegrand = figure; scatter(xlat, y, 10)
    title(sprintf('%s_n-%d_Tx-%s', ...
      fName, npts, ptransform), 'interpreter','none')
    saveas(hFigIntegrand, sprintf('%s_n-%d_Tx-%s_rFun-%1.2f.png', ...
      fName, npts, ptransform, rfun))
  end
  
%   if 0
%     lnaMLE_opt = fminsearch(@(lna) ...
%       ObjectiveFunction(exp(lna),r,xlat,(ftilde)), ...
%       0,optimset('TolX',1e-2));
%     thetaOpt = exp(lnaMLE_opt)
%   else
%     if 0
%       thetaOpt = 1;
%       ln_rOpt = fminsearch(@(lnr) ...
%         ObjectiveFunction(thetaOpt,1+exp(lnr),xlat,(ftilde)), ...
%         0,optimset('TolX',1e-2));
%       rOpt = 1 + exp(ln_rOpt)
%       r = rOpt;
%     else %search for optimal kernel parameters
      objfun = @(lnParams) ...
        ObjectiveFunction(exp(lnParams(1)),1+exp(lnParams(2)),xlat,(ftilde));
      %% Plot the objective function
      lnthetarange = (-2:0.2:2);
      lnorderrange = (-1:0.1:1);
      [lnthth,lnordord] = meshgrid(lnthetarange,lnorderrange);
      objobj = lnthth;
      for ii =  1:size(lnthth,1)
         for jj = 1:size(lnthth,2)
            objobj(ii,jj) = objfun([lnthth(ii,jj); lnordord(ii,jj)]);
         end
      end
      figure
      s = surf(lnthth,lnordord,objobj);
      set(s,'EdgeColor','none','facecolor','interp')
      set(gca,'xtick',log([0.2 0.4 1 3 7]),'xticklabel',{'0.2','0.4','1','3','7'}, ...
         'ytick',log([1.4 1.6 2 2.6 3.7]-1), 'yticklabel',{'1.4', '1.6','2','2.6','3.7'})
      xlabel('\(\theta\)')
      ylabel('\(r\)')
      
      [~,which] = min(objobj,[],'all','linear');
      [whichrow,whichcol] = ind2sub(size(lnthth),which);
      lnthOptAppx = lnthth(whichrow,whichcol);
      thetaOptAppx = exp(lnthOptAppx)
      lnordOptAppx = lnordord(whichrow,whichcol);
      orderOptAppx = 1+exp(lnordOptAppx)

      
      %% Optimize the objective function
      lnParamsOpt = fminsearch(objfun,[lnthOptAppx;lnordOptAppx],optimset('TolX',1e-3));
      thetaOpt = exp(lnParamsOpt(1));
      rOpt = 1 + exp(lnParamsOpt(2));
      hold on 
      scatter3(lnParamsOpt(1),lnParamsOpt(2),objfun(lnParamsOpt)*1.002,1000,MATLABYellow,'.')%     end
%   end
  
  % lambda1 = kernel(r, xlat_, thetaOpt);
  vlambda = kernel2(thetaOpt, rOpt, xlat);
  s = sqrt(sum(abs(ftilde(2:end).^2)./vlambda(2:end))/(npts^2));
  vlambda = (s^2)*vlambda;

  % apply transform
  % $\vZ = \frac 1n \mV \mLambda^{-\frac 12} \mV^H(\vf - m \vone)$
  % ifft also includes 1/n division
  vz = ifft(ftilde./sqrt(vlambda));
  vz_real = real(vz);  % vz must be real as intended by the transformation 

  % create_plots('normplot')
  create_plots('qqplot', vz_real, fName, npts, ptransform, r, rOpt, theta, thetaOpt)
  
  % Shapiro-Wilk test
  %   [H, pValue, W] = swtest(w_ftilde);
  %   Hval='true';
  %   if H==true
  %     Hval='false';
  %   end
  %   fprintf('Shapiro-Wilk test: Normal=%s, pValue=%1.3f, W=%1.3f\n', Hval, pValue, W);
  fprintf('r = %7.5f, rOpt = %7.5f, theta = %7.5f, thetaOpt = %7.5f\n', ...
     r,rOpt,theta,thetaOpt);
  

end

end

function create_plots(type, vz_real, fName, npts, ptransform, r, rOpt, theta, thetaOpt)
hFigNormplot = figure();
set(hFigNormplot,'defaultaxesfontsize',16, ...
  'defaulttextfontsize',16, ... %make font larger
  'defaultLineLineWidth',0.75, 'defaultLineMarkerSize',8)
if strcmp(type, 'normplot')
  normplot(vz_real)
else
  qqplot(vz_real); hold on
  plot([-3 3], [-3 3],'-')
end

title(sprintf('%s r=%1.2f rOpt=%1.2f, theta=%1.2f, thetaOpt=%1.2f', ...
       fName, r, rOpt, theta, thetaOpt))
saveas(hFigNormplot, sprintf('%s_%s_n-%d_Tx-%s_rOpt-%1.3f.png', ...
       type, fName, npts, ptransform, r))
end

% gaussian random function
function fval = f_rand(xpts, rfun, a, b, c, seed)
dim = size(xpts,2);
rng(seed) % initialize random number generator for reproducability
N1 = 2^floor(16/dim);
Nall = N1^dim;
kvec(dim,Nall) = 0; %initialize kved
kvec(1,1:N1) = 0:N1-1; %first dimension
Nd = N1;
for d = 2:dim
   Ndm1 = Nd;
   Nd = Nd*N1;
   kvec(1:d,1:Nd) = [repmat(kvec(1:d-1,1:Ndm1),1,N1); ...
      reshape(repmat(0:N1-1,Ndm1,1),1,Nd)];
end
f_c = a*randn(1,Nall-1);
f_s = a*randn(1,Nall-1);
f_0 = c + b*randn(1, 1);
kbar = prod(max(kvec(:,2:Nall),1),1);
argx = (2*pi*xpts) * kvec(:,2:Nall);
f_c_ = (f_c./kbar.^(rfun)).*cos(argx);
f_s_ = (f_s./kbar.^(rfun)).*sin(argx);
fval = f_0 + sum(f_c_+ f_s_,2) ;
% figure; plot(fval, '.')
end

function Lambda = kernel(theta, r, xun)
constMult =  -(-1)^(r/2)*((2*pi)^r)/(2*factorial(r));
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

function vlambda = kernel2(theta, r, xun)
n = size(xun, 1);
m = (1:(-1 + n/2))';
tilde_g_h1 = m.^(-r);
tilde_g = [0; tilde_g_h1; 0; tilde_g_h1(end:-1:1)];
g = fft(tilde_g);
temp_ = (theta/2)*g(1 + xun*n);
C1 = prod(1 + temp_, 2);
% matlab's builtin fft is much faster and accurate
% eigenvalues must be real : Symmetric pos definite Kernel
vlambda = real(fft(C1));
end

% function vlambda = kernel2(r, xun, theta)
% constMult = 1/2;
% kernelFunc = @(x) truncated_series_kernel(x,r);
%         
% temp_ = (theta*constMult)*kernelFunc(xun);
% C1 = prod(1 + temp_, 2);
% 
% % matlab's builtin fft is much faster and accurate
% % eigenvalues must be real : Symmetric pos definite Kernel
% vlambda = real(fft(C1));
% end
% function g = truncated_series(N, r)
% tilde_g_0 = 0;
% m = 1:(-1 + N/2);
% tilde_g_h1 = N./abs(m).^(r);
% m = (N/2):(-1 + N);
% tilde_g_h2 = N./abs(N-m).^(r);
% tilde_g = [tilde_g_0 tilde_g_h1 tilde_g_h2];
% g = ifft(tilde_g)';
% end
% 
% function c = truncated_series_kernel(x,r)
% n = size(x, 1);
% g = truncated_series(n,r);
% c = g(1 + x*n);
% end

function [loss,Lambda,RKHSnorm] = ObjectiveFunction(theta,order,xun,ftilde)
tol = 100*eps;
n = length(ftilde);
% [Lambda, Lambda_ring] = kernel(xun,obj.kernType,a,obj.order,...
%   obj.avoidCancelError);
arbMean = true;
[Lambda] = kernel2(theta ,order, xun);

% compute RKHSnorm
%temp = abs(ftilde(Lambda~=0).^2)./(Lambda(Lambda~=0)) ;
temp = abs(ftilde(Lambda>tol).^2)./(Lambda(Lambda>tol)) ;

% compute loss: MLE
if arbMean==true
  RKHSnorm = sum(temp(2:end))/n;
  temp_1 = sum(temp(2:end));
else
  RKHSnorm = sum(temp)/n;
  temp_1 = sum(temp);
end

% ignore all zero eigenvalues
loss1 = sum(log(Lambda(Lambda>tol)))/n;
loss2 = log(temp_1);
loss = (loss1 + loss2);
if imag(loss) ~= 0
   keyboard
end
%fprintf('L1 %1.3f L2 %1.3f L %1.3f r %1.3e theta %1.3e\n', loss1, loss2, loss, order, theta)

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

