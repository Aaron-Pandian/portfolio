function u_backwards_Euler = Backwards_Euler_3(A, u0, f, n, T, dt)

Nt = round(T/dt);
u_backwards_Euler = zeros(n, Nt+1);
u_backwards_Euler(:,1) = u0;

a = zeros(1, n);
b = zeros(1, n);
c = zeros(1, n);
d = zeros(n, n);

for i = 1:n
    d(i,i) = 1;
end

f2 = d - dt*A;

for i = 1:n
    if i > 1
        a(i) = f2(i, i-1);
    end
    
    b(i) = f2(i,i);
    
    if i < n
        c(i) = f2(i, i+1);
    end
end

for i = 1:Nt
    r = u_backwards_Euler(:,i) + dt*f;
    u_backwards_Euler(:,i+1) = Tridiag(a,b,c,r);
end
end
