clear all
clear clc
%% b = 10
lambda = 1;
figure(1);
x = linspace(-20,20);
y = linspace(-.25,.25);
[X,Y] = meshgrid(x,y);
v = 10;
b = 10;
theta1 = atan2(Y,(X+b));
theta2 = atan2(Y,(X-b));
Z = v*Y + (lambda/(2*pi))*(theta1-theta2);
contour(X,Y,Z);
[C,h] = contour(X,Y,Z);
xlabel("x");
ylabel("y");
title("b = 10")
clabel(C,h,0);
%% b = 1
lambda = 1;
figure(2);
x = linspace(-2,2);
y = linspace(-.25,.25);
[X,Y] = meshgrid(x,y);
v = 10;
b = 1;
theta1 = atan2(Y,(X+b));
theta2 = atan2(Y,(X-b));
Z = v*Y + (lambda/(2*pi))*(theta1-theta2);
contour(X,Y,Z);
[C,h] = contour(X,Y,Z);
xlabel("x");
ylabel("y");
title("b = 1")
clabel(C,h,0);
%% b = .01
lambda = 1;
figure(3);
x = linspace(-.2,.2);
y = linspace(-.2,.2);
[X,Y] = meshgrid(x,y);
v = 10;
b = .01;
theta1 = atan2(Y,(X+b));
theta2 = atan2(Y,(X-b));
Z = v*Y + (lambda/(2*pi))*(theta1-theta2);
contour(X,Y,Z);
[C,h] = contour(X,Y,Z);
xlabel("x");
ylabel("y");
title("b = .01")
clabel(C,h,0);
%% b = .001
lambda = 1;
figure(4);
x = linspace(-.03,.03);
y = linspace(-.03,.03);
[X,Y] = meshgrid(x,y);
v = 10;
b = .001;
theta1 = atan2(Y,(X+b));
theta2 = atan2(Y,(X-b));
Z = v*Y + (lambda/(2*pi))*(theta1-theta2);
contour(X,Y,Z);
[C,h] = contour(X,Y,Z);
xlabel("x");
ylabel("y");
title("b = .001");
clabel(C,h,0);
%% b = .0001
lambda = 1;
figure(5);
x = linspace(-.015,.015);
y = linspace(-.015,.015);
[X,Y] = meshgrid(x,y);
v = 10;
b = .0001;
theta1 = atan2(Y,(X+b));
theta2 = atan2(Y,(X-b));
Z = v*Y + (lambda/(2*pi))*(theta1-theta2);
contour(X,Y,Z);
[C,h] = contour(X,Y,Z);
xlabel("x");
ylabel("y");
title("b = .0001");
clabel(C,h,0);
%% t/c vs b 
b = [.0001,.001,.01,1,10];
tc1 = (lambda/(pi*v*.0001))*((pi/2)-atan(.005/.0001));
tc2 = (lambda/(pi*v*.001))*((pi/2)-atan(.017/.001));
tc3 = (lambda/(pi*v*.01))*((pi/2)-atan(.3/.01));
tc4 = (lambda/(pi*v*1))*((pi/2)-atan(.45/1));
tc5 = (lambda/(pi*v*10))*((pi/2)-atan(.5/10));
tc = [tc1,tc2,tc3,tc4,tc5];
figure(6);
plot(b,tc);
xlabel("b");
ylabel("t/c");
%{ 
% Why the Rakine Oval becomes a cylinder as b approaches zero: 
% At high values of b, the oval is streched thin with a chord length much larger than the thickness. 
% But at extremely low values of b, the length of the chord and the thickness are nearly equal. 
% Thus, as b decreases, the chord length and thickness even out, thus the value of t/c approaches 1. 
% A t/c value of 1 represents a circle or clylinder.
%}
%% Lambda Source = .1, Lambda Sink = 1
b = 1;
v = 10;
lambdasink = 1;
lambdasource = .1;
figure(7);
x = linspace(-5,5);
y = linspace(-.5,.5);
[X,Y] = meshgrid(x,y);
theta1 = atan2(Y,(X+b));
theta2 = atan2(Y,(X-b));
Z = v*Y + (lambdasource/(2*pi))*theta1-(lambdasink/(2*pi))*theta2;
contour(X,Y,Z);
[C,h] = contour(X,Y,Z);
xlabel("x");
ylabel("y");
title("Sink Stronger");
clabel(C,h,0);
%% Lambda Source = 1, Lambda Sink = .1
b = 1;
v = 10;
lambdasink = .1;
lambdasource = 1;
figure(8);
x = linspace(-5,5);
y = linspace(-.5,.5);
[X,Y] = meshgrid(x,y);
theta1 = atan2(Y,(X+b));
theta2 = atan2(Y,(X-b));
Z = v*Y + (lambdasource/(2*pi))*theta1-(lambdasink/(2*pi))*theta2;
contour(X,Y,Z);
[C,h] = contour(X,Y,Z);
xlabel("x");
ylabel("y");
title("Source Stronger");
clabel(C,h,0);







