function J = CostFunction3(x, params, Fx)
% Write function that computes cost function J
% 
% x is a vector of optimization parameters
% params is a struct that contains all other parameters
% F is a function that takes x as input and returns solution of tumor model
% get g and f from tumor model

J = 0;
delta=x(1:params.Ntreat);

[g, f, t] = Fx(x);

deltaT = t(2) - t(1);

for i = 1:params.Ntreat
    J = J + delta(i)^2;
end

J = J + params.a*g(end);

for i = 1:length(g)-1
    J = J + params.b*((g(i) + g(i+1))/2 *deltaT);
end

for i = 1:length(f)-1
    J = J + params.c*((f(i) + f(i+1))/2 *deltaT);
end

