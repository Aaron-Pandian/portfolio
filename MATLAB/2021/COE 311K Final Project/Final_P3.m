%% Regression
L2_Step = 0.05;
L2_Values = zeros(16,1);
h_function_values = zeros(16,1);

for i = 0:15
    L2_Values(i+1) = 0.02 + i*L2_Step;
    [h_function_values(i+1),A] = Derrivation_Function(L2_Values(i+1));
end

regression = polyfit(L2_Values, h_function_values, 13);
regression_function = @(x) polyval(regression, x) - 50;
L2_regression = Bisection(0.07, 0.1, regression_function);

Final_Table = [L2_Values, h_function_values];

approximated_h = Derrivation_Function(L2_regression);
error = (50 - approximated_h)/50;

%% Optimization
f = @(L2) Derrivation_Function(L2);
J = @(L2) (f(L2)-50)^2;
% Equation given in problem set

tolx = 1.e-10;
tolfun = 1.e-10;
initial_guess = 0.07;
iterations = 1000;
minimum = 0.05;
maximum = 0.2;

L2_optimized = fmincon(J,initial_guess,[],[],[],[],minimum,maximum,[],optimset('TolX',tolx,'TolFun',tolfun,'MaxIter', iterations,'Display','iter-detailed'));
