function [H_Values, A_New]= Derrivation_Function(L2)
%% Parameters and Initial Conditions

L12 = 0.1;
L1 = 0.7;
L = L1 + L12 + L2;
T = 1;
Nt = 1000;
dt = T/Nt;
Nx = 1000;
dx = L/Nx;
ht0 = 100; 
k1 = 10;
k2 = 0.1; 
alpha = 200;
beta = 2;
theta = 20;
c = 0.001;

%% Initilizing A Matrices 

A = zeros(Nx, Nx);
A(1, 1) = -2/(dx^2);
A(1, 2) = 1/(dx^2);

for row = 2:Nx-1
    A(row, row-1) = 1/(dx^2);
    A(row, row) = -2/(dx^2);
    A(row, row+1) = 1/(dx^2);
end

index = Nx-1;
A(index+1, index) = 1/(dx^2);
A(index+1, index+1) = -1/(dx^2);

k = @(x) K_Func(x, k1, k2, L1, L12);

for i = 1:Nx
    row = round(i);
    A(row,:) = k(row*dx) .* A(row,:);
end

A2 = zeros(Nx, Nx);
% Matrix for second partial derivative

for row = 1:Nx-1
    A2(row,row) = (-1/dx);
    A2(row,row+1) = (1/dx);
end

A_New = A - c*A2;

%% Initilizing Qext
Lx = 0:dx:L;
qext = @(x) alpha.*exp(-beta.*x).*sin(theta.*x);
Qext = qext((Lx)');
Qext = Qext(2:Nx+1);
Qext(1,1) = Qext(1,1) + k(dx)*(ht0)/(dx^2);

%% Solving using Backwards Euler
n = Nx;
u0 = 0;
u_Backwards_Euler = Backwards_Euler_3(A_New, u0, Qext, n, T, dt);
u_Backwards_Euler_Table = [ht0; u_Backwards_Euler(:,Nt+1)];
H_Values = u_Backwards_Euler_Table(Nx+1);

end