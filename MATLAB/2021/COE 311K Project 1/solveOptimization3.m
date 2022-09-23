function[g3, f3]=solveOptimization3()

clear all; clc; clf;
% Set parameters 
% e.g.
delta = [.001,.001,.001,.001];
tau = [9.7955,19.5171,30.6384,40.6257];
x=[delta, tau];

sigma = 2; a = .0005; b = .00005; c = .001; Ntreat = 4;

T = 50; k = 4; %error for k4=0.2167

Lp = .1; La = .005; Lk = 4; Ld = .05;

g0 = .01; f0 = .001;

% add other parameters
% create Matlab struct to collect all parameters in one object
params.T = T; params.g0 = g0; params.f0 = f0; params.tau = tau; 
params.sigma = sigma; params.a = a; params.b = b; params.c = c; 
params.Ntreat = Ntreat; params.La = La; params.Lp = Lp;
params.Lk = Lk; params.Ld = Ld; params.k = k;params.x=x;

% keep adding other parameters to 'params'
% you can access values from 'params' using params.T, params.g0, ...
% Optimization setup
% set initial guess for optimization parameters

% define function that computes cost function J
% using struct 'params' you don't need to pass lengthy list of parameters
Fx = @(x) TumModel3(x, params);  % tumor model
Jx = @(x) CostFunction3(x, params, Fx); % cost function


minTau=zeros(1,4); 
maxTau=zeros(1,4);
for i=1:4
    minTau(i)= 5+10*(i-1);
    maxTau(i)=15+10*(i-1);
end

% set constraints (lower and upper bounds)
xmin = [0,0,0,0,minTau]; % vector of lower bound for all elements of x vector
xmax = [0.01,0.01,0.01,0.01,maxTau]; % vector of upper bound for all elements of x vector

% fmincon parameters (you may leave this unchanged)
tolx = 1.e-9;
tolfun = 1.e-9;
maxiter = 400;
% Run optimization

tic;

[xopt, Jval, ~, ~, ~, ~, ~] = fmincon(Jx, x, [], [], [], [], ...
                                xmin, xmax,[], ...
                                optimset('TolX',tolx, ...
                                'TolFun', tolfun, ...
                                'MaxIter', maxiter, ...
                                'Display','iter-detailed'));

toc;

% display results
disp('initial guess for delta')
disp(delta)
disp('optimal parameters for delta')
disp(xopt(1:4))
disp('initial guess for tau')
disp(tau)
disp('optimal parameters for tau')
disp(xopt(5:8))

% Plotting
[g3, f3, ~]=Fx(xopt);
