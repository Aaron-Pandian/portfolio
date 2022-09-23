function [y] = backward_euler2(A, h0, T, dt, q, n)
h= h0;
I = eye(n);

for i = 1:(T/dt)
    J = (I - (dt*A));
    b = h(:,i)+dt*q;
    h(:,i+1) = J\b;
end

y = h;
end