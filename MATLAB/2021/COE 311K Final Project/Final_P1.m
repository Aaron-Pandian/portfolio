n = 4;
T = 1;
dt = 10^-5;
f = [0,0,0,0];
u0 = [1;0;0;0];
A = [-2,1,0,0; 1,-2,1,0;0,1,-2,1;0,1,-1,0];

%% Part 1 and 3 Forward Euler Method
disp('Part 1 & 3, Forward Euler Method u(T)')
forward = forward_euler(A,u0,T,dt,f);
disp(forward(:,4))

%% Part 2 and 3 Backward Euler Method
disp('Part 2 & 3, Backward Euler Method u(T)')
backward = backward_euler(A,u0,T,dt,f,n);
disp(backward)

%% Part 4
T_new = 1000;
dt_new = T_new/10^6;

disp('Part 4, Forward Euler u(T) for large T')
forward_new = forward_euler(A,u0,T_new,dt_new,f);
disp(forward_new(:,4))

disp('Part 4, Backward Euler u(T) for large T')
backward_new = backward_euler(A,u0,T_new,dt_new,f,n);
disp(backward_new)

%% Part 5
T2 = 10;
dt2_values = zeros(1,10)';

for k = 1:10
    
    dt2 = T/(2^k);
    dt2_values(k) = dt2;
    u = backward_euler(A,u0,T2,dt2,f,n);
    maxU(k) = max(u);
    
end

k_values = (1:1:10)';
table = [k_values, dt2_values, (maxU)'];
disp('Part 5, Max(u(T)) the maximum of four-element vector u(T) for different k from 1 to 10')
disp('Column 1 are k Values, Column 2 are the delta t values, and Column 3 are the max u(T) values')
disp(table)