function [y] = forward_euler2(dhdt, u0, T, dt)
u = u0;

for i = 2:(T/dt)
    u(:,i) = u(:,i-1) + dt*dhdt(u(:,i-1));
end

y = u;
end