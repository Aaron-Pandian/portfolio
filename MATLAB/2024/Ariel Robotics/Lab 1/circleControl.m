quadParamsScript;
constantsScript;

r = [0.5 0 0.5]'; % not editable
v = [0 -1 0]'; % must always be in the Y direction only 
radius = 1;

m = quadParams.m;
g = constants.g;

% set euler matrix 
forceRatio = -(m*(v(2)^2))/(m*g*radius);
theta = atan(forceRatio);
phi = 0;
psi = 0;
e = [phi, theta, psi]; % e(1) = phi, e(2) = theta, e(3) = psi

% find wB matrix
yawRate = (2*pi)/((2*pi*radius)/v(2));
angularRates = [0 0 yawRate]'; % negative yaw for clockwise rotation around origin
conRotOmegaMat = [cos(phi)*cos(theta) 0 cos(phi)*cos(theta);
                  sin(phi)*sin(theta) cos(phi) -cos(theta)*sin(phi);
                  -sin(theta) 0 cos(theta)];
wB = conRotOmegaMat*angularRates;

% find NB
J = quadParams.Jq;
NB = crossProductEquivalent(wB)*J*wB;

% find w12 using w34, where wb = w34 is less than wa = w12 
% 607 for speed 2
wb = 585;
kF = quadParams.kF(1);
appliedForce = -(m*(v(2)^2))/(radius*sin(theta)); % put negative becuase centripital in negative x dir
Fz = appliedForce*cos(theta);
Fx = appliedForce*sin(theta);
Fzwa = 2*kF*(wb^2)*cos(theta);
Fzwb = Fz - Fzwa;
wa = sqrt(Fzwb/(2*kF*cos(theta)));

fprintf('wB = %f\n',wB')
fprintf('e = %f\n',e')
fprintf('w12 = %f, w34 = %f\n',wa,wb)
Fxapproximation = 2*kF*(wa^2)*sin(theta) + 2*kF*(wb^2)*sin(theta);
fprintf('Fx = %f, Fxp = %f\n',Fx,Fxapproximation)









