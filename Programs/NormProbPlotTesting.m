%% Normality testing
clearvars, close all
set(0,'defaultLineMarkerSize',10)

rng default;  % For reproducibility
n = 16;

%%
x = (0:(n-1))'/n;
V = exp((2*pi/n)*sqrt(-1)*(0:(n-1))'*(0:(n-1)));
VH = transpose(conj(V));
max(V*transpose(conj(V))-n*eye(n))
C = 1 - mod(x - x',1) + (mod(x - x',1)).^2;
Lambda =  VH*C(:,1);
Cmhalf = V*diag(Lambda.^(-1/2))*VH;
max(max(abs(imag(Cmhalf))))
return

%% Some plots
mu = 10;
sigma = 2;
X = normrnd(mu,sigma,n,1); %normal, but not standard normal data

figure; normplot(X)

Z = (X - mean(X))/std(X); %standardized normal data

disp(['+/- one standard deviation = [' ...
   num2str(normcdf(-1)) ',' num2str(normcdf(1)) ']']) 

figure; normplot(Z)

figure; qqplot(Z)

figure; qqplot(X)

Y = norminv((1/n : 1/n : 1 - 1/n)'); %ideal standard normal quantiles

figure; qqplot(Z,Y); hold on; plot([-4 4],[-4 4],'g')

figure; qqplot(X,Y); hold on; plot([-4 4],[-4 4],'g')

%% Complex data
A = randn(n,n) + sqrt(-1)*randn(n,n);
[V,~,~] = svd(A);
max(abs(V*transpose(conj(V))) - eye(n))

Zcomp = V*Z;
mean(Zcomp)
std(Zcomp)
Zreal = real(Zcomp);
Zimag = imag(Zcomp);
Zall = [Zreal; Zimag];
Yall = norminv((1/(2*n) : 1/(2*n) : 1 - 1/(2*n))');
figure; qqplot(Zall,Yall); hold on; plot([-4 4],[-4 4],'g')



