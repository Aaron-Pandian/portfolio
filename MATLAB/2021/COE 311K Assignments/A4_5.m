a = 0;
b = 3;
syms x;
syms y;
f1 = x*exp(1.5*x);
intf_exact = double(int(f1,x,a,b));
disp('The exact integration of f is:')
disp(intf_exact);

% Gauss Quadrature 
f = @(x) x.*exp(1.5.*x); 
x_loc = [-1/sqrt(3) 1/sqrt(3)];
x_loc2 = [-sqrt(3/5) 0 sqrt(3/5)];
w = [1,1];  
w2 = [(5/9),(8/9),(5/9)];
g = @(y)(b-a).*f(a+(b-a).*(y+1)./2)./2;
value = 0;
for i = 1:2
    value = value + w(i)*g(x_loc(i));
end
disp('The value for the two-point method is:')
disp(value);
error1 = 100*(intf_exact-value)/intf_exact;
disp('The error for the two-point method is:')
disp(error1)

value2 = 0;
for i=1:3
    value2=value2 + w2(i)*g(x_loc2(i));
end
disp('The value for the three-point method is:')
disp(value2);
error2 = 100*(intf_exact-value2)/intf_exact;
disp('The error for the three-point method is:')
disp(error2)
