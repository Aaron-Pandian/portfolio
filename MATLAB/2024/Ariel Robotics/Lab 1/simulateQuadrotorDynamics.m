function [P] = simulateQuadrotorDynamics(S)
% simulateQuadrotorDynamics : Simulates the dynamics of a quadrotor aircraft.
%
%
% INPUTS
%
% S ---------- Structure with the following elements:
%
%          tVec = Nx1 vector of uniformly-sampled time offsets from the
%                 initial time, in seconds, with tVec(1) = 0.
%
%  oversampFact = Oversampling factor. Let dtIn = tVec(2) - tVec(1). Then the
%                 output sample interval will be dtOut =
%                 dtIn/oversampFact. Must satisfy oversampFact >= 1.   
%
%      omegaMat = (N-1)x4 matrix of rotor speed inputs.  omegaMat(k,j) >= 0 is
%                 the constant (zero-order-hold) rotor speed setpoint for the
%                 jth rotor over the interval from tVec(k) to tVec(k+1).
%
%        state0 = State of the quad at tVec(1) = 0, expressed as a structure
%                 with the following elements:
%                   
%                   r = 3x1 position in the world frame, in meters
% 
%                   e = 3x1 vector of Euler angles, in radians, indicating the
%                       attitude
%
%                   v = 3x1 velocity with respect to the world frame and
%                       expressed in the world frame, in meters per second.
%                 
%              omegaB = 3x1 angular rate vector expressed in the body frame,
%                       in radians per second.
%
%       distMat = (N-1)x3 matrix of disturbance forces acting on the quad's
%                 center of mass, expressed in Newtons in the world frame.
%                 distMat(k,:)' is the constant (zero-order-hold) 3x1
%                 disturbance vector acting on the quad from tVec(k) to
%                 tVec(k+1).
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
% P ---------- Structure with the following elements:
%
%          tVec = Mx1 vector of output sample time points, in seconds, where
%                 P.tVec(1) = S.tVec(1), P.tVec(M) = S.tVec(N), and M =
%                 (N-1)*oversampFact + 1.
%                  
%  
%         state = State of the quad at times in tVec, expressed as a structure
%                 with the following elements:
%                   
%                rMat = Mx3 matrix composed such that rMat(k,:)' is the 3x1
%                       position at tVec(k) in the world frame, in meters.
% 
%                eMat = Mx3 matrix composed such that eMat(k,:)' is the 3x1
%                       vector of Euler angles at tVec(k), in radians,
%                       indicating the attitude.
%
%                vMat = Mx3 matrix composed such that vMat(k,:)' is the 3x1
%                       velocity at tVec(k) with respect to the world frame
%                       and expressed in the world frame, in meters per
%                       second.
%                 
%           omegaBMat = Mx3 matrix composed such that omegaBMat(k,:)' is the
%                       3x1 angular rate vector expressed in the body frame in
%                       radians, that applies at tVec(k).
%
%
%+------------------------------------------------------------------------------+
% References:
%
%
% Author: Aaron Pandian
%+==============================================================================+  

% Initializing vectors
tVec = S.tVec;
oversampFact = S.oversampFact;
omegaMat = S.omegaMat;
state0 = S.state0;
distMat = S.distMat;
r = state0.r;
e = state0.e;
v = state0.v;
omegaB = state0.omegaB;

% counting all bases by the subsequent, seemingly nonsencical formula
% for creating X0 
rI = r;
vI = v;
RBI = euler2dcm(e); % transformation
X0 = [rI(1,1), rI(2,1), rI(3,1), vI(1,1), vI(2,1), vI(3,1), RBI(1,1), RBI(2,1), RBI(3,1), RBI(1,2), RBI(2,2), RBI(3,2), RBI(1,3), RBI(2,3), RBI(3,3), omegaB(1,1), omegaB(2,1), omegaB(3,1)]';

% Oversampling causes the ODE solver to produce output at a finer time
% resolution than dtIn. oversampFact is the oversampling factor.
N = length(tVec);
dtIn = tVec(2) - tVec(1);
dtOut = dtIn/oversampFact;

% Create empty storage vectors
tVecOut = [];
XMat = [];

% Set initial state 
Xk = X0;

% Starting the ode45 method
for k = 1:N-1
  
    % Create instances of the  quad motor dynamics function parameters
    distVec = distMat(k,:)';
    omegaVec = [omegaMat(k,1) omegaMat(k,2) omegaMat(k,3) omegaMat(k,4)]';
    
    % Build the time vector for kth segment.  We oversample by a factor
    % oversampFact relative to the coarse timing of each segment.
    tspan = [tVec(k):dtOut:tVec(k+1)]';
    
    % Run ODE solver for time segment
    [tVeck,XMatk] = ode45(@(t,X)quadOdeFunction(t,X,omegaVec,distVec,S), tspan, Xk);
    
    % Add the data from the kth segment to your storage vector
    tVecOut = [tVecOut; tVeck(1:end-1)];
    XMat = [XMat; XMatk(1:end-1,:)];
    
    % Prepare for the next iteration
    Xk = XMatk(end,:)';
end

% Store the final state of the final segment
XMat = [XMat; XMatk(end,:)];
tVecOut = [tVecOut; tVeck(end,:)];

% Updating N
N = length(tVecOut);

% Creating P Structure
P.tVec = tVecOut; % size N x 1?
rMatrix = zeros(3, N); 
eMatrix = zeros(3, N);
vMatrix = zeros(3, N);
omegaBmatrix = zeros(3, N);

for i = 1:N % iterations on time so matrix should look like 18 x X
    Xi = XMat';
    rMatrix(1:end,i) = [Xi(1,i), Xi(2,i), Xi(3,i)]'; 
    vMatrix(1:end,i) = [Xi(4,i), Xi(5,i), Xi(6,i)]';
    RBIk = [Xi(7,i) -Xi(8,i) Xi(9,i); 
            Xi(10,i) Xi(11,i) Xi(12,i); 
            Xi(13,i) Xi(14,i) Xi(15,i)]';
    ek = dcm2euler(RBIk); % transform back to euler angles for output
    ek(3) = -ek(3); % accounting for negative error
    eMatrix(1:end,i) = ek;
    omegaBmatrix(1:end,i) = [Xi(16,i) Xi(17,i) Xi(18,i)]';
end

outputState.rMat = rMatrix';
outputState.eMat = eMatrix';
outputState.vMat = vMatrix';
outputState.omegaBMat = omegaBmatrix';
P.state = outputState;





  

