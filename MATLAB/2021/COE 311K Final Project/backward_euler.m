function [y] = backward_euler(A, u0, T, dt, f, n)
u = u0;
I = eye(n);

for i = 1:(T/dt)
    J = (I - (dt*A));
    u = J\u;
end

y = u;
end