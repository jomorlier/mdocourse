%% Introduction to Kriging (Gaussian Process)
%% 				Prof. Joseph MORLIER
% Figure 1 illustrates a typical example of a prediction problem: given some 
% noisy observations of a dependent variable at certain values of the independent 
% variable , what is our best estimate of the dependent variable at a new value 
% x* ? 

%clear all; close all;
printoptions;


%  Inputs
X = [-1.5 -1 -.75 -.4 -.25 0]'; n = size(X,1);
% Outputs
y = .55*[-3 -2 -.6 .4 1 1.6]'; 
%
figure;
plot(X,y,'k.')
xlabel('X'); ylabel('y');

%% 
% What if we want to estimate the model at a new point  x* defined as 18 
% points in [-1.5,0]. Please check also x*= 0.2 (extrapolation) !!
% 
% The regression problem is defined as follows:
% 
% Let $$\mathbf{x}_i \in {R}^{6}$ $ be an input vector and $$\mathbf{y}_i 
% \in {R}^{6}$$ be its corresponding target. The set of $$M$$ inputs are arranged 
% into a matrix $$\mathbf{X} = [\mathbf{x}_1, \dots, \mathbf{x}_M]^\top$$ and 
% their corresponding targets are stored in a matrix $$\mathbf{Y} = [\mathbf{y}_1 
% - \mathbf{\bar{y}}, \dots, \mathbf{y}_M-\mathbf{\bar{y}}]^\top$, with $\mathbf{\bar{y}}$ 
% $being the mean target value in $$\mathbf{Y}$$.
% 
% 
% 
% 

%xstar= 0.2; %extrapolation
xstar=-1.5:0.1: 0.2; %18 points in [-1.5,0] interpolation
N = length(xstar);
%computing covariance in matrix form
[covXXInd1 covXXInd2] = meshgrid(X,X)

%computing Kstar
[covXXsInd1 covXXsInd2] = meshgrid(X,xstar)
%% 
% 		

%computing Kstarstar
[covXsXsInd1 covXsXsInd2] = meshgrid(xstar,xstar)
%% 
% 				
% 
% hey hey add something missing, we may need parameters to fit "a model" 
% of our covariance matrix. Let's try with a Standard Exponential (SE) Kernel.
% 
% $$k(x,x') =\sigma_f^2\exp\left(-\frac{(x-x')^2}{l^ 2}\right)$$
% 
% 

l = 1; % first hyperparameter
sig_f = sqrt(3); % second

covXsXs = sig_f * exp(-(covXsXsInd1-covXsXsInd2).^2 ./ l.^2)
covXX = sig_f * exp(-(covXXInd1-covXXInd2).^2 ./ l.^2)

covXXs = sig_f * exp(-(covXXsInd1-covXXsInd2).^2 ./ l.^2)
figure;
imagesc(covXXs); % Plotting the Gram Matrix
%% 
% I wish to train a GPR model $$\mathcal{G} = \lbrace \mathbf{X}, \mathbf{Y}, 
% \theta \rbrace$$ using the squared exponential function ($$\theta$$ must be 
% chosen):
% 
% 
% 
% $$k(\mathbf{x}_i, \mathbf{x}_j) = \sigma_f^2 \text{exp} \left( - \frac{1}{l^2}(\mathbf{x}_i 
% - \mathbf{x}_j)^2\right) + \sigma_n^2\delta_{ij}$$, where $$\delta_{ij}$$ equals 
% $$1$$ if $$i = j$$ and $$0$$ otherwise. 

sig_n=0.3
 %sig_n=0 %instead

covXX_noisy = covXX + sig_n^2 *eye(size(covXX))
figure;
imagesc(covXX_noisy); % Plotting the Gram Matrix

posterior_mean = covXXs/covXX_noisy * y 
figure;plot(xstar,posterior_mean); hold on;
plot(X,y,'k.'); xlabel('X'); ylabel('y');
%% 
% Compute the 95% confidence interval
% 
% 

posterior_cov = covXsXs - covXXs/covXX_noisy * covXXs'
%% 
% 		

% Plotting mean, variance and sample of the posterior % Plotting the variance of Prior
f = [posterior_mean+2*sqrt(diag(posterior_cov));flipdim(posterior_mean -2*sqrt(diag(posterior_cov)),1) ];
%%bounds = [posterior_mean+1.96*sqrt(posterior_cov) posterior_mean-1.96*sqrt(posterior_cov)]

xf=[xstar'; flipdim(xstar',1)];
figure; fill(xf, f, [7 7 7]/8);hold on; plot(xstar,posterior_mean); hold on;
plot(X,y,'k.'); xlabel('X'); ylabel('y');
%% 
% 	

% Tip: add a jitter term to the gram matrix so that matrix inversion is numerically stable
jitter = 10^(-6);
% The cholesky decomposition
L = chol(posterior_cov + eye(N)*jitter);
% Plotting the randomly drawn sample function
randomFunction = posterior_mean + L'*randn(N,5); 
figure; fill(xf, f, [7 7 7]/8);hold on; plot(xstar , randomFunction);
%% 
%  but is there a way to find the $\theta_{optimal}$ ? oh $\theta = [l, 
% \sigma_f, \sigma_n]$
% 
% The hyperparameters are $\theta = [l, \sigma_f, \sigma_n]$ with $$\sigma_n$$ 
% being the assumed noise level in the training data and $$l$ $ is the length-scale  
% (of oscillations) and $$\sigma_f$$ the amplitude.
% 
% To train the model, I need to minimise the negative log marginal likelihood 
% with respect to the hyperparameters:
% 
% $$$-\text{log}\, p(\mathbf{Y} \mid \mathbf{X}, \theta) = \frac{1}{2} \text{tr}(\mathbf{Y}^\top\mathbf{K}^{-1}\mathbf{Y}) 
% + \frac{1}{2}\text{log}\mid\mathbf{K}\mid + \,c,$$$
% 
% where c is a constant and the matrix $$\mathbf{K}$ $is a function of the 
% hyperparameters (see equation k(xi,xj) = ...).
% 
% Based on the demo given on the GPML website, Let's try to implement this 
% using the GPML Matlab code is below.
% 
%