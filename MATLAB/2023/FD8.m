Cd = .025;
Cl = .349;
Fnstatic = 216000;
Kf = Fnstatic;
Ka0 = 1;
Ka1 = 3.281*(10^-5);
Ka2 = 10.764*(10^-9);
Km0 = 1;
Km1 = -1.07;
Km2 = .56;
U = .03;
W0 = 790100;
Wfdot = 9;
S = 125;
g = 9.81;
p = 1.225;
L = .5*p*S*Cl;
D = .5*p*S*Cd;
h = 156.1;
f1 = Ka0 + Ka1*h + Ka2*(h^2);
f2 = 1;
P0 = 101325;
Papt = 99534;
Pratio = Papt/P0;
F = Pratio*Kf*f1*f2;

%% Part A EOM for Takeoff Roll

% Euler Method Numerical Integration
% https://www.mathworks.com/matlabcentral/answers/278300-matlab-code-help-on-euler-s-method

h = (.05);  % step size
t = (0):h:(24);  % the range of x
V = zeros(size(t));  % allocate the result y
V(1) = (0);  % the initial y value
At = zeros(size(t));  % allocate the result y
At(1) = (0);  % the initial y value
Dt = zeros(size(t));  % allocate the result y
Dt(1) = (0);  % the initial y value
n = numel(V);  % the number of y values
% The loop to solve the DE
for i=1:n-1
    A = (g/(W0-(Wfdot*t(i)))) * (F*(Ka0 + Ka1*V(i) + Ka2*((V(i))^2)) - (D*((V(i))^2)) - U*((W0-(Wfdot*t(i))) - L*((V(i))^2)));
    At(i) = A;
    if i > 1
        Dt(i) = Dt(i-1) + h * V(i);
    end
    V(i+1) = V(i) + h * A;
end

figure(1)
% Acceleration
plot(t,At)
title("Acceleration V. Time");
xlabel('Time (s)') 
ylabel('Acceleration (m^2/s)') 
% Velocity
figure(2)
plot(t,V)
title("Velocity V. Time");
xlabel('Time (s)') 
ylabel('Velocity (m/s)')
% Displacement 
figure(3)
plot(t,Dt)
title("Displacement V. Time");
xlabel('Time (s)') 
ylabel('Distance (m)')

%% Part B  

% Tweaked equation for constant weight 
Wfdot = 0;

h = (.05);  % step size
t = (0):h:(24);  % the range of x
V = zeros(size(t));  % allocate the result y
V(1) = (0);  % the initial y value
Dt = zeros(size(t));  % allocate the result y
Dt(1) = (0);  % the initial y value
n = numel(V);  % the number of y values
% The loop to solve the DE
for i=1:n-1
    A = (g/(W0-(Wfdot*h))) * (F*(Ka0 + Ka1*V(i) + Ka2*((V(i))^2)) - (D*((V(i))^2)) - U*((W0-(Wfdot*h)) - L*((V(i))^2)));
    if i > 1
        Dt(i) = Dt(i-1) + h * V(i);
    end
    V(i+1) = V(i) + h * A;
end

% Velocity
figure(4)
plot(t,V)
% Displacement 
figure(5)
plot(t,Dt)

%% Part C 

% Tweaked equation for constant weight and acceleration constant
Wfdot = 0;

h = (.05);  % step size
t = (0):h:(24);  % the range of x
V = zeros(size(t));  % allocate the result y
V(1) = (0);  % the initial y value
Dt = zeros(size(t));  % allocate the result y
Dt(1) = (0);  % the initial y value
n = numel(V);  % the number of y values
% The loop to solve the DE
A = (g/(W0-(Wfdot*h))) * (F*(Ka0 + Ka1*V(i) + Ka2*((V(i))^2)) - (D*((V(i))^2)) - U*((W0-(Wfdot*h)) - L*((V(i))^2)));
for i=1:n-1
    if i > 1
        Dt(i) = Dt(i-1) + h * V(i);
    end
    V(i+1) = V(i) + h * A;
end

% Velocity
figure(6)
plot(t,V)
% Displacement 
figure(7)
plot(t,Dt)







