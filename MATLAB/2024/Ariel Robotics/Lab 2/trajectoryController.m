function [Fk,zIstark] = trajectoryController(R,S,P)
% trajectoryController : Controls quadcopter toward a reference trajectory.
%
%
% INPUTS
%
% R ---------- Structure with the following elements:
%
%       rIstark = 3x1 vector of desired CM position at time tk in the I frame,
%                 in meters.
%
%       vIstark = 3x1 vector of desired CM velocity at time tk with respect to
%                 the I frame and expressed in the I frame, in meters/sec.
%
%       aIstark = 3x1 vector of desired CM acceleration at time tk with
%                 respect to the I frame and expressed in the I frame, in
%                 meters/sec^2.
%
% S ---------- Structure with the following elements:
%
%        statek = State of the quad at tk, expressed as a structure with the
%                 following elements:
%                   
%                  rI = 3x1 position in the I frame, in meters
% 
%                 RBI = 3x3 direction cosine matrix indicating the
%                       attitude
%
%                  vI = 3x1 velocity with respect to the I frame and
%                       expressed in the I frame, in meters per second.
%                 
%              omegaB = 3x1 angular rate vector expressed in the body frame,
%                       in radians per second.
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
% Fk --------- Commanded total thrust at time tk, in Newtons.
%
% zIstark ---- Desired 3x1 body z axis expressed in I frame at time tk.    
%                  
%+------------------------------------------------------------------------------+
% References:
%
%
% Author:  
%+==============================================================================+  

m = P.quadParams.m;
g = P.constants.g;

rIstark = R.rIstark;
vIstark = R.vIstark;
aIstark = R.aIstark;

rI = S.statek.rI;
RBI = S.statek.RBI;
vI = S.statek.vI;

er = rIstark - rI;
er_dot = vIstark - vI;

K = [4 0 0; 0 4 0; 0 0 4]; % Parameters to tune
Kd = [1.5 0 0; 0 1.5 0; 0 0 1.5]; % Both 3x3

FIstark = K*er + Kd*er_dot + [0; 0; m*g] + m*aIstark;
zIstark = FIstark./norm(FIstark);

e3 = [0 0 1]';
zI = (RBI')*e3;
Fk = (FIstark')*zI;







  

