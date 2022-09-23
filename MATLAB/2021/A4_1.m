clear all;
T = 10;
dt = 0.001;
m = 1;
c = 0.05;
k = 0.75;
u0 = 0;
udot0 = 1;

%% Part 1
t = 0:dt:T;

% Velocity-Verlet method
f = @(u,v) (-c/m)*v - (k/m)*u;
u_verlet = V_Verlet(T, dt, u0, udot0, f);
hold on;
plot(t, u_verlet,'r','DisplayName','Verlet')

% Forward Euler method
steps = round(T/dt);
u = zeros(1, steps+1);
u2 = zeros(1, steps+1);
v = zeros(1, steps+1);
x = zeros(1, steps+1);
u(:,1) = u0;
u2(:,1) = u0;
v(:,1) = udot0;
f = @(u,v) -(c./m).*v - (k./m).*u;
f1 = @(u2,v) -(c./m).*v - (k./m).*u2;
f2 = @(u2,x) -(c./m).*x - (k./m).*u2;
for i = 1:steps
    u(:,i+1) = u(:,i) + dt*v(:,i);
    v(:,i+1) = v(:,i) + dt*f(u(i),v(i));
end
uT = u(:, steps);
vT = v(:, steps);

plot(t, u,'b','DisplayName','Forward Euler')
legend()

%% Part 2
k = 0.75;
delta_k = 0.01;
N = 50;
sigma = 0.1;
k_values = zeros(N+1,1);
uT_table = zeros(N+1, 10);

k_values(1) = k;
f = @(u,v) (-c/m)*v - (k_values(1)/m)*u;
for i = 1:10
    x_rand = sigma*randn(1,1);
    u_rand = V_Verlet(T, dt, x_rand, udot0, f);
    uT = u_rand(length(u_rand));
    uT_table(1,i) = uT;
end

for i = 2:N+1
    k_values(i) = k_values(i-1) + delta_k;
    f = @(u,v) (-c./m).*v - (k_values(i)./m).*u;
    for j = 1:10
        x_rand = sigma*randn(1,1);
        u_rand = V_Verlet(T, dt, x_rand, udot0, f);
        uT = u_rand(length(u_rand));
        uT_table(i,j) = uT;
    end
end

answer = [k_values, uT_table];
%disp(answer)

figure;
plot(k_values, uT_table, '.', 'Color', 'black')
hold on;

%% Part 4
% Chose Regression and Polynomial 
x4 = 0.5:.001:1.5;

u_k = polyfit(k_values*ones(1,10), uT_table, 2);
f4 = polyval(u_k, x4);
plot(x4, f4, 'red')

%{
u_k = polyfit(k_values*ones(1,10), uT_table, 13);
f4 = polyval(u_k, x4);
plot(x4, f4, 'green')
%}

%{
u_k = polyfit(k_values*ones(1,10), uT_table, 9);
f4 = polyval(u_k, x4);
plot(x4, f4, 'yellow')
%}

%% Part 5
k_vals = [0.65, 0.71, 0.83, 0.96, 1.02,1.09,1.17,1.26,1.34]';
uT_old = zeros(length(k_vals), 1);
uT_new = zeros(length(k_vals), 1);
for i = 1:length(k_vals)
    uT_old(i) = polyval(u_k, k_vals(i));
    f = @(u,v) (-c./m).*v - (k_vals(i)./m).*u;
    u_new = V_Verlet(T, dt, u0, udot0, f);
    uT_new(i) = u_new(length(u_new));
end
error = (uT_old - uT_new)/100;
answer2 = [k_vals, uT_old, uT_new, error];
disp(answer2);
