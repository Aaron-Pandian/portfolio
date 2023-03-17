disp("NEW");
% Shock undergoes 3 stages, with 3 turn angles, find property relations
% Given ma and theta2 = theta_tot - theta1
mach = 2.6;
gamma = 1.4;
theta_tot = 20;
theta1deg = 9;
theta1 = theta1deg*pi/180;
theta2 = (theta_tot-theta1deg)*pi/180;
theta3 = (theta2);

%% Oblique Shock Calculations - Step 1 
% Change m0norm to m1norm, and change theta to next number ie theta2..
m0norm = mach;
theta=theta1; 

% Find Wave Angle - θ-β-M Equation reference linked below 
mu=asin(1/m0norm);                   
c=tan(mu)^2;
n = 0;
a=((gamma-1)/2+(gamma+1)*c/2)*tan(theta);
b=((gamma+1)/2+(gamma+3)*c/2)*tan(theta);
d=sqrt(4*(1-3*a*b)^3/((27*a^2*c+9*a*b-2)^2)-1);
beta=atan((b+9*a*c)/(2*(1-3*a*b))-(d*(27*a^2*c+9*a*b-2))/(6*a*(1-3*a*b))*tan(n*pi/3+1/3*atan(1/d)))*180/pi;
% Rudd, L., and Lewis, M. J., "Comparison of Shock Calculation Methods", AIAA Journal of Aircraft, Vol. 35, No. 4, July-August, 1998, pp. 647-649.
beta = (beta)*(pi/180);

% Find minitial oblique
m0ob = m0norm*sin(beta);
% Calculate mfinal oblique 
m1ob = ((1+((gamma-1)/2)*(m0ob^2))/((gamma*(m0ob^2))-((gamma-1)/2)))^.5; 
% Calculate mfinal normal 
m1norm = (m1ob)/(sin(beta-theta)); 
% Calculate static pressure ratios p1/pa, p2/p1, p3/p2, p3/pa
staticpratio = 1+(((2*gamma)/(gamma+1))*((m0ob^2)-1));
% Calculate density ratio
staticdratio = ((gamma+1)*(m0ob^2))/(2+((gamma-1)*(m0ob^2)));
% Calculate temperature ratio
statictempratio = staticpratio*(1/staticdratio);
% Calculate stagnation pressure ratios p2o/p1o, p3o/p2o, p3o/pao
stagpratio = ((staticdratio)^(gamma/(gamma-1)))/((staticpratio)^(1/(gamma-1)));

% Output necessary value (for homework): 
x = [m1norm, stagpratio, staticpratio];
display(x);

%% Calculations - Step 2 
m0norm = m1norm;
theta=theta2; 

% Find Wave Angle - Equation reference linked below 
mu=asin(1/m0norm);                   
c=tan(mu)^2;
n = 0;
a=((gamma-1)/2+(gamma+1)*c/2)*tan(theta);
b=((gamma+1)/2+(gamma+3)*c/2)*tan(theta);
d=sqrt(4*(1-3*a*b)^3/((27*a^2*c+9*a*b-2)^2)-1);
beta=atan((b+9*a*c)/(2*(1-3*a*b))-(d*(27*a^2*c+9*a*b-2))/(6*a*(1-3*a*b))*tan(n*pi/3+1/3*atan(1/d)))*180/pi;
% Rudd, L., and Lewis, M. J., "Comparison of Shock Calculation Methods", AIAA Journal of Aircraft, Vol. 35, No. 4, July-August, 1998, pp. 647-649.
beta = (beta)*(pi/180);

% Find minitial oblique
m0ob = m0norm*sin(beta);
% Calculate mfinal oblique 
m1ob = ((1+((gamma-1)/2)*(m0ob^2))/((gamma*(m0ob^2))-((gamma-1)/2)))^.5; 
% Calculate mfinal normal 
m1norm = (m1ob)/(sin(beta-theta)); 
% Calculate static pressure ratios p1/pa, p2/p1, p3/p2, p3/pa
staticpratio2 = 1+(((2*gamma)/(gamma+1))*((m0ob^2)-1));
% Calculate density ratio
staticdratio2 = ((gamma+1)*(m0ob^2))/(2+((gamma-1)*(m0ob^2)));
% Calculate stagnation pressure ratios p2o/p1o, p3o/p2o, p3o/pao
stagpratio2 = ((staticdratio2)^(gamma/(gamma-1)))/((staticpratio2)^(1/(gamma-1)));

% Output necessary value: 
x = [m1norm, stagpratio2, staticpratio2];
display(x);

%% Calculations - Step 3 (Normal Shock)
m0norm = m1norm;
theta=theta3; 

% Calculate mfinal normal 
m1norm = ((1+((gamma-1)/2)*(m0norm^2))/((gamma*(m0norm^2))-((gamma-1)/2)))^.5; 
% Calculate static pressure ratios p1/pa, p2/p1, p3/p2, p3/pa
staticpratio3 = 1+(((2*gamma)/(gamma+1))*((m0norm^2)-1));
% Calculate density ratio
staticdratio3 = ((gamma+1)*(m0norm^2))/(2+((gamma-1)*(m0norm^2)));
% Calculate stagnation pressure ratios p2o/p1o, p3o/p2o, p3o/pao
stagpratio3 = ((staticdratio3)^(gamma/(gamma-1)))/((staticpratio3)^(1/(gamma-1)));

% Output necessary value: 
x = [m1norm, stagpratio3, staticpratio3];
display(x);

%% Calculations - Step 4 Ambient 
stagpratio4 = stagpratio*stagpratio2*stagpratio3;
staticpratio4 = staticpratio*staticpratio2*staticpratio3;

% Output necessary value: 
x = [mach, stagpratio4, staticpratio4];
display(x);


