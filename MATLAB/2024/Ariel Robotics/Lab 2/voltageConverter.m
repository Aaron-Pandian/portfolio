function [eak] = voltageConverter(Fk,NBk,P)
% voltageConverter : Generates output voltages appropriate for desired
%                    torque and thrust.
%
%
% INPUTS
%
% Fk --------- Commanded total thrust at time tk, in Newtons.
%
% NBk -------- Commanded 3x1 torque expressed in the body frame at time tk, in
%              N-m.
%
% P ---------- Structure with the following elements:
%
%    quadParams = Structure containing all relevant parameters for the
%                 quad, as defined in quadParamsScript.m 
%
%     constants = Structure containing constants used in simulation and
%                 control, as defined in constantsScript.m 
%
%
% OUTPUTS
%
% eak -------- Commanded 4x1 voltage vector to be applied at time tk, in
%              volts. eak(i) is the voltage for the ith motor.
%
%+------------------------------------------------------------------------------+
% References:
%
%
% Author:  
%+==============================================================================+  

locations = P.quadParams.rotor_loc;
kF = P.quadParams.kF(1);
kN = P.quadParams.kN(1);
kT = kN/kF;
Cm = P.quadParams.cm(1);
eaMax = P.quadParams.eamax;
beta = 0.9; % constant
alpha = 1; % can reduce

% Finding max force value
wmax = Cm*eaMax;
Fmax = kF*(wmax^2);
FmaxTotal = Fmax*4*beta;
extraVecOne = [FmaxTotal, Fk];

G = [1, 1, 1, 1;
     locations(2,1), locations(2,2), locations(2,3), locations(2,4);
     -locations(1,1), -locations(1,2), -locations(1,3), -locations(1,4);
     -kT, kT, -kT, kT];

% Expressed by min force and torque vector
extraVec = [min(extraVecOne), alpha*NBk(1), alpha*NBk(2), alpha*NBk(3)]';
Fvec = (G^(-1))*extraVec;

% changing alpha to satisfy Fi <= Fmax 
while max(all(Fvec > Fmax)) % if any element is over Fmax return a 1 to indicate true, if all are false, max is a 0 or false. 
    alpha = alpha - 0.05;
    extraVec = [min(extraVecOne), alpha*NBk(1), alpha*NBk(2), alpha*NBk(3)]';
    Fvec = (G^(-1))*extraVec; % recalculate Fvec to check again, if passes, this is new value
end

% Check if motor force values are below zero and fix
for F = 1:4
    if Fvec(F) < 0
        Fvec(F) = 0;
    end
end

omegaVec = ((1/kF)*Fvec).^(0.5);
eak = (1/Cm)*omegaVec;







  

