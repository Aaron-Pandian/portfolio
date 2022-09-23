clear all;
% Problem Constants and vectors
L12 = 0.1; a = 200; b = 2; theta = 20; T = 1; L = 1; Nt = 1000; dt = T/Nt; 
Nx = 1000; dx = L/Nx; k1 = 10; k2 = 0.1; L1 = 0.7; L2 = 0.2; 
xvec = 0:dx:L; tvec = 0:dt:T; n = length(xvec)-1;
Qext =  a*exp(-b*xvec).*sin(theta*xvec);
Htimederiv = [-100/dx.^2; zeros(length(tvec)-2,1)]';
H = [100 zeros(1,length(xvec)-1); zeros(length(tvec)-1,length(xvec))];
H(:,1) = 100; hz0 = H(1:end-1,1); h0 = zeros(Nx,1);
% Creates k vector
for i = 1:length(xvec)
    if xvec(i) <= L1
        k(i) = k1;
    elseif xvec(i) >= L1 + L12
        k(i) = k2;
    else
        k(i) = k1 + (k2 - k1) *(xvec(i) - L1)/L12;
    end
end

% making second derivative matrix A

%creates 1 -2 1 diagnal matrix

e = ones(n,1);
A = spdiags([e -2*e e],-1:1,n,n);
A(end,end) = -1;
for i = 1:Nx
    A(i,:) = A(i,:)*k(i)/(dx^2);
end
A = full(A);

%Creates q matrix
for i = 1:(length(tvec)-1)
    for j = 1:(length(xvec))
        q(i,j) = a * exp(-b*xvec(j))*sin(theta*xvec(j));
    end
end

q(:,1) = q(:,1) + 10000 * (k(i)/(dx^2));
q = [q(:,1),q(:,3:end)];
q = q(1,:)';

%Finds h using backward Euler
h = backward_euler2(A,h0,T,dt,q,n);

%Plots h at k values
ktest = [201,401,601,801,1001];
size(h(:,201))
size(xvec)

plot(xvec, [100;h(:,ktest(1))])
hold on
plot(xvec, [100;h(:,ktest(2))])
hold on
plot(xvec, [100;h(:,ktest(3))])
hold on
plot(xvec, [100;h(:,ktest(4))])
hold on
plot(xvec, [100;h(:,ktest(5))])
hold on
legend("201","401","601","801","1001")

hold off


%Part3
dhdt = @(H) A*H+q;

T = .001;
Nt = 20000;
dt = T/Nt;
hfe = forward_euler2(dhdt, h0, T, dt);

% Finds best Nt value for Forward Euler
% Nt = (5000:1000:50000);
% 
% for i = 1:length(Nt)
%     dt = T/Nt(i);
%     hfe = forward_euler2(dhdt, h0, T, dt);
%     out(i) = hfe(end,end);
% end
%
% 5 X 10^(-8) is the largest dt we can have that still outputs a stable dataset. 