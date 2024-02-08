function [Xdot] = quadOdeFunction(t,X,omegaVec,distVec,P)
% quadOdeFunction : Ordinary differential equation function that models
%                   quadrotor dynamics.  For use with one of Matlab's ODE
%                   solvers (e.g., ode45).
%
%
% INPUTS
%
% t ---------- Scalar time input, as required by Matlab's ODE function
%              format.
%
% X ---------- Nx-by-1 quad state, arranged as 
%
%              X = [rI',vI',RBI(1,1),RBI(2,1),...,RBI(2,3),RBI(3,3),omegaB']'
%
%              rI = 3x1 position vector in I in meters
%              vI = 3x1 velocity vector wrt I and in I, in meters/sec
%             RBI = 3x3 attitude matrix from I to B frame
%          omegaB = 3x1 angular rate vector of body wrt I, expressed in B
%                   in rad/sec
%
% omegaVec --- 4x1 vector of rotor angular rates, in rad/sec.  omegaVec(i)
%              is the constant rotor speed setpoint for the ith rotor.
%
%  distVec --- 3x1 vector of constant disturbance forces acting on the quad's
%              center of mass, expressed in Newtons in I.
%
% P ---------- Structure with the following elements:
%
%    quadParams = Structure containing all relevant parameters for the
%                 quad, as defined in quadParamsScript.m 
%
%     constants = Structure containing constants used in simulation and
%                 control, as defined in constantsScript.m 
%
% OUTPUTS
%
% Xdot ------- Nx-by-1 time derivative of the input vector X
%
%+------------------------------------------------------------------------------+
% References:
%
%
% Author: Aaron Pandian
%+==============================================================================+

% Initalizing matricies and variables
m = P.quadParams.m;
g = P.constants.g;
KF = P.quadParams.kF;
KN = P.quadParams.kN;
J = P.quadParams.Jq; 
ri = P.quadParams.rotor_loc; % 3 X 4
dI = distVec; % 3 x 1
rI = [X(1,1) X(2,1) X(3,1)]'; % 3 x 1, not needed for subsequent calculations
vI = [X(4,1) X(5,1) X(6,1)]'; % 3 x 1

RBI = [X(7,1) X(10,1) X(13,1); % 3 x 3
       X(8,1) X(11,1) X(14,1); 
       X(9,1) X(12,1) X(15,1)]; 

wB = [X(16,1) X(17,1) X(18,1)]'; % 3 x 1
wi = omegaVec; % 4 x 1

% Finding rI_dot
rI_dot = vI;

% Finding vI_dot
FzI = [0 0 -m*g]';
F1 = [0 0 KF(1,1)*((wi(1,1))^2)]'; % Zeroed to the z-axis due to formula
F2 = [0 0 KF(2,1)*((wi(2,1))^2)]';
F3 = [0 0 KF(3,1)*((wi(3,1))^2)]';
F4 = [0 0 KF(4,1)*((wi(4,1))^2)]';
rotor_force_term = F1 + F2 + F3 + F4;
vI_dot = (1/m) * (FzI +  (transpose(RBI)*rotor_force_term) + dI);

% Finding w_dot
r1 = [ri(1,1) ri(2,1) ri(3,1)]'; % 3 X 1
r2 = [ri(1,2) ri(2,2) ri(3,2)]';
r3 = [ri(1,3) ri(2,3) ri(3,3)]';
r4 = [ri(1,4) ri(2,4) ri(3,4)]';
N1 = -1*[0 0 KN(1)*((wi(1,1))^2)]'; % The 1 and 3 propellor have negative values from formula
N2 =  1*[0 0 KN(2)*((wi(2,1))^2)]'; % The 2 and 4 propellor have positive values from formula
N3 = -1*[0 0 KN(3)*((wi(3,1))^2)]';
N4 =  1*[0 0 KN(4)*((wi(4,1))^2)]';

NB1 = N1 + crossProductEquivalent(r1)*F1;
NB2 = N2 + crossProductEquivalent(r2)*F2;
NB3 = N3 + crossProductEquivalent(r3)*F3;
NB4 = N4 + crossProductEquivalent(r4)*F4;
NB = NB1 + NB2 + NB3 + NB4;

w_dot = (inv(J)) * (NB - crossProductEquivalent(wB)*J*wB); % Replace inv() method

% Finding RBI_dot 
RBI_dot = -(crossProductEquivalent(wB))*RBI;

% Initializing output
Xdot = [rI_dot(1,1) rI_dot(2,1) rI_dot(3,1) vI_dot(1,1) vI_dot(2,1) vI_dot(3,1) RBI_dot(1,1) RBI_dot(2,1) RBI_dot(3,1) RBI_dot(1,2) RBI_dot(2,2) RBI_dot(3,2) RBI_dot(1,3) RBI_dot(2,3) RBI_dot(3,3) w_dot(1,1) w_dot(2,1) w_dot(3,1)]';





