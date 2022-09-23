x = linspace(-1,1,11);
y = 1./(1+25*x.^2);
x2 = linspace(-1,1);
f = polyfit(x,y,10);
y10 = polyval(f,x2);
yr = 1./(1+25*x2.^2);
plot(x,y,'--',x2,y10,x2,yr,'g')
n = 11;
k = 1:11;
f = @(x) 1./(1+25*x.^2); 
xc = cos(((2*k-1)*pi)/22);
yc_nodes = f(xc);

hold on;
plot(xc,yc_nodes,'b')

% Polynomial Derrivation 
f = polyfit(xc,yc_nodes,n-1);
plot(x,polyval(f,x),'r'); 
legend('Chebyshev Nodes','Lagrunge Function','Interpolating Polynomial')