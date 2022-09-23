function [y] = forward_euler(A, u0, T, dt, f)
u = u0;

for i = 1:(T/dt)
    u = u + dt*(A*u+f);
end

y = u;
end