function[g1, f1, g2, f2, t]=solveOptimization()
clear all; clc; clf;
% Set parameters 
% e.g.
delta = [.001,.001,.001,.001];
tau = [10,20,30,40];

sigma = 2; a = .0005; b = .00005; c = .001; Ntreat = 4;

T = 50; k = 4; %error for k4=0.2167

Lp = .1; La = .005; Lk = 4; Ld = .05;

g0 = .01; f0 = .001;

% add other parameters
% create Matlab struct to collect all parameters in one object
params.T = T; params.g0 = g0; params.f0 = f0; params.tau = tau; 
params.sigma = sigma; params.a = a; params.b = b; params.c = c; 
params.Ntreat = Ntreat; params.La = La; params.Lp = Lp;
params.Lk = Lk; params.Ld = Ld; params.k = k;

% keep adding other parameters to 'params'
% you can access values from 'params' using params.T, params.g0, ...
% Optimization setup
% set initial guess for optimization parameters

% define function that computes cost function J
% using struct 'params' you don't need to pass lengthy list of parameters
Fx = @(x) TumModel(x, params);  % tumor model
Jx = @(x) CostFunction(x, params, Fx); % cost function


error = drugError(delta, params);
% set constraints (lower and upper bounds)
xmin = [0,0,0,0];  % vector of lower bound for all elements of x vector
xmax = [0.01,0.01,0.01,0.01];  % vector of upper bound for all elements of x vector

% fmincon parameters (you may leave this unchanged)
tolx = 1.e-9;
tolfun = 1.e-9;
maxiter = 400;
% Run optimization

tic;

[xopt, Jval, ~, ~, ~, ~, ~] = fmincon(Jx, delta, [], [], [], [], ...
                                xmin, xmax,[], ...
                                optimset('TolX',tolx, ...
                                'TolFun', tolfun, ...
                                'MaxIter', maxiter, ...
                                'Display','iter-detailed'));

toc;

% display results
disp('initial guess')
disp(delta)
disp('optimal parameters')
disp(xopt)
% Plotting
 [g1,f1,t]=Fx(delta);
 [g2,f2,~]=Fx(xopt);

subplot(2,1,1)
plot(t, g1)
hold on
plot(t, g2)
legend('gi','go')

subplot(2,1,2)
plot(t, f1)
hold on 
plot(t, f2)
legend('fi','fo')
