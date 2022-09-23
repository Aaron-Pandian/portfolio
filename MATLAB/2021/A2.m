%% Problem 1 Part 1 
syms x
syms qext(x)
qext(x) = 12*x^2 + cos(5*x) + 100*x*sin(10*x);
disp('external heat function')
disp(qext)
syms Q1(x)
Q1(x) = int(qext, 0, x);
disp('Q1 function')
disp(Q1)
syms Q2(x)
Q2(x) = int(Q1, 0, x);
disp('Q2 function')
disp(Q2)
disp('1.1) value of Q2 function at x = 1')
constant = (Q2(1));
disp(constant)

%% Problem 1 Part 2 
x = 0:.1:1;
qext = @(x) 12*x.^2 + cos(5*x) + 100*x.*sin(10*x);
T = @(x) 100*x + 200*x*(36/25 - sin(10) - (2*cos(5)^2)/5 - cos(5)/25) - (x.^4 - x.*sin(10*x) -2*cos(5*x).*cos(5*x)/5 - cos(5*x)/25 + 11/25)*200;
plot(x, qext(x), 'r+', 'DisplayName', 'External Heat')
set(gcf, 'Name', 'Question 1.2')
hold on;
plot(x, T(x), 'r', 'DisplayName', 'Temperature')
hold on;
y1 = 0*x + 80;
y2 = 0*x + 40;
plot(x, y1, 'g', 'DisplayName', 'Line y = 80')
hold on;
plot(x, y2, 'k', 'DisplayName', 'Line y = 40')
legend()

%% Problem 2 Part 1 
x = linspace(0,1,6);
t = @T;
for i = x
    if t(i) <= 80 && t(i+.2) >= 80 
        x1 = i; 
        x2 = i+.2; 
    end
end
answer = ['2.1) The interval is ',num2str(x1),' - ',num2str(x2),'.'];
disp(answer)

%% Problem 2 Part 2
x_lower = x1;
x_upper = x2;

x_mid = (x_lower + x_upper) / 2;

while abs(Tnew(x_mid))>.00001
   if Tnew(x_mid) * Tnew(x_upper) < 0
       x_lower = x_mid;
   else
       x_upper = x_mid;
   end
   x_mid = (x_lower + x_upper) / 2;
end

fprintf('2.2) The root approximation is %g\n',x_mid)

%% Problem 3 Part 1 
Q1 = @(x) sin(5*x)/5 + sin(10*x)+10*x*(2*sin(5*x)^2 - 1)+ 4*x^3;
Q2 = @(x) x^4 - x*sin(10*x) -2*cos(5*x)*cos(5*x)/5 - cos(5*x)/25 + 11/25 ;
dt = @(x) 100 + 200*(Q2(1)) - 200*(Q1(x));
iter = 0;
xi = .15;
tol = .00001;
maxiter = 100;

while(1)
   xiold = xi;
   xi = xi - (Tlast(xi)/dt(xi));
   iter = iter + 1;
  if xi ~= 0, err = abs((xi-xiold)/xi)*100; end
  if err <= tol || iter >= maxiter, break, end
end
fprintf('3.1) The root approximation is %g\n',xi)    
    
    


%% Part 3 Part 2 
x1 = .15;
h = .0001;
x2 = x1 - h;
tolerance = .00001;
f1 = Tlast(x1);
dx = inf;
iter = 0;

while abs(dx) > tolerance && iter <= 100
    iter = iter + 1;
    f2 = Tlast(x2);
    dx = (x2 - x1)*f2/(f2-f1);
    x1 = x2;
    f1 = f2;
    x2 = x2 - dx;
end
thisanswer = ['3.2) The root approximation is ',num2str(x2),'.'];
disp(thisanswer)

%% Question 4 Part 1
x_lower = .6;
x_upper = .9;
x_mid = (x_lower + x_upper) / 2;
Q1 = @(x) sin(5*x)/5 + sin(10*x)+10*x*(2*sin(5*x).^2 - 1)+ 4*x.^3;
Q2 = @(x) x.^4 - x.*sin(10*x) -2*cos(5*x).*cos(5*x)/5 - cos(5*x)/25 + 11/25 ;
f = @(x) 100 + 200*(Q2(1)) - 200*(Q1(x));

while abs(f(x_mid))>.00001
   if f(x_mid) * f(x_upper) < 0
       x_lower = x_mid;
   else
       x_upper = x_mid;
   end
   x_mid = (x_lower + x_upper) / 2;
end

max = T(x_mid);
fprintf('4.1) The max temp approximation is at %g\n',x_mid)
fprintf('The max temp approximation is %g degrees celcius\n',max)

%% Question 4 Part 2
x_0 = .6;
x_1 = .9;
x = 0:.01:1;
Q1 = @(x) sin(5*x)/5 + sin(10*x)+10*x*(2*sin(5*x).^2 - 1)+ 4*x.^3;
Q2 = @(x) x.^4 - x.*sin(10*x) -2*cos(5*x).*cos(5*x)/5 - cos(5*x)/25 + 11/25 ;
T =@(x)100*x+ 200*x.*(Q2(1))-200.*(Q2(x));
f = @(x) 100 + 200*(Q2(1)) - 200*(Q1(x));
err = abs(x_1 - x_0);
root = 0;
while err>.000001
   x2=(x_0*f(x_1)-x_1*f(x_0))/(f(x_1)-f(x_0));
    x_0=x_1;
    x_1=x2;
    err=abs(x_1-x_0);
    root=x2;
end

max = T(root);
fprintf('4.2) The maximum tempererature is at %g\n',root)
fprintf('So the maximum temperature would be %g degrees celcius\n',max)





