function [NBk] = attitudeController(R,S,P)
% attitudeController : Controls quadcopter toward a reference attitude
%
%
% INPUTS
%
% R ---------- Structure with the following elements:
%
%       zIstark = 3x1 desired body z-axis direction at time tk, expressed as a
%                 unit vector in the I frame.
%
%       xIstark = 3x1 desired body x-axis direction, expressed as a
%                 unit vector in the I frame.
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
% NBk -------- Commanded 3x1 torque expressed in the body frame at time tk, in
%              N-m.
%
%+------------------------------------------------------------------------------+
% References:
%
%
% Author:  
%+==============================================================================+  

J = P.quadParams.Jq;

RBI = S.statek.RBI;
wB = S.statek.omegaB;
% Small angle assumption
eE_dot = wB;

zIstark = R.zIstark;
xIstark = R.xIstark;

K = [1 0 0; 0 1 0; 0 0 1]; % Parameters to tune
Kd = [0.3 0 0; 0 0.3 0; 0 0 0.3];

b = cross(zIstark, xIstark)/(norm(cross(zIstark, xIstark))); % 3x1
a = cross(b,zIstark); % 3x1
RBIstark = [a, b, zIstark]'; % 3x3
RE = RBIstark*(RBI'); % Double check if element by element multiplication is needed here
eE = [RE(2,3)-RE(3,2), RE(3,1)-RE(1,3), RE(1,2)-RE(2,1)]'; % 3x1

NBk = K*eE - Kd*eE_dot + crossProductEquivalent(wB)*J*wB;

  

