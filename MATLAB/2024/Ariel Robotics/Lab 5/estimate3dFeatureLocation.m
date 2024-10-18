function [rXIHat,Px] = estimate3dFeatureLocation(M,P)
% estimate3dFeatureLocation : Estimate the 3D coordinates of a feature point
%                             seen by two or more cameras with known pose.
%
%
% INPUTS
%
% M ---------- Structure with the following elements:
%
%       rxArray = 1xN cell array of measured positions of the feature point
%                 projection on the camera's image plane, in pixels.
%                 rxArray{i} is the 2x1 vector of coordinates of the feature
%                 point as measured by the ith camera.  To ensure the
%                 estimation problem is observable, N must satisfy N >= 2 and
%                 at least two cameras must be non-colinear.
%
%      RCIArray = 1xN cell array of I-to-camera-frame attitude matrices.
%                 RCIArray{i} is the 3x3 attitude matrix corresponding to the
%                 measurement rxArray{i}.
%
%       rcArray = 1xN cell array of camera center positions.  rcArray{i} is
%                 the 3x1 position of the camera center corresponding to the
%                 measurement rxArray{i}, expressed in the I frame in meters.
%
% P ---------- Structure with the following elements:
%
%  sensorParams = Structure containing all relevant parameters for the quad's
%                 sensors, as defined in sensorParamsScript.m
%
% OUTPUTS
%
%
% rXIHat -------- 3x1 estimated location of the feature point expressed in I
%                 in meters.
%
% Px ------------ 3x3 error covariance matrix for the estimate rxIHat.
%
%+------------------------------------------------------------------------------+
% References:
%
%
% Author:  
%+==============================================================================+ 

% Define parameters
Rc_value = P.sensorParams.Rc(1,1);
Rc = [Rc_value, Rc_value];
K = P.sensorParams.K;
pixelSize = P.sensorParams.pixelSize;

% rxArray{i} is the 2x1 vector of the feature point as measured by the ith camera.
rx = M.rxArray;

% RCIArray{i} is the 3x3 attitude matrix 
RCI = M.RCIArray;

% rcArray{i} is the 3x1 position of the camera center in I frame
tI = M.rcArray;

H = [];
RcPrime = [];
% Begin for loop through cameras to create H and R
for i = 1:length(rx)

    % Rotation from Oc to Oi
    RCIk = RCI{i};
    
    % Translation from Oi to Oc in C
    tIk = -tI{i};
    tCk = RCIk*tIk;

    % Projection matrix for ith camera
    Pk = K*[RCIk tCk];

    % Initialize projection point in meters
    xk = pixelSize*rx{i};

    % Update H
    Pkr1 = Pk(1,:);
    Pkr2 = Pk(2,:);
    Pkr3 = Pk(3,:);
    Hk = [xk(1)*Pkr3 - Pkr1; xk(2)*Pkr3 - Pkr2];
    H = [H; Hk];

    % Update RcPrime array
    RcPrime = [RcPrime, Rc]; % Creating array of length 2N
end

% Create covariance matrix R
RcPrime = diag(RcPrime);
R = (pixelSize^2)*RcPrime;

% Create Hr and z matrices
Hr = H(:,1:end-1);
z = -1*H(:,end);

% Solve for error covarience matrix
Px = (Hr'*(R^(-1))*Hr)^(-1);

% Rearrange to solve for Xr, location of feature point 
rXIHat = Px*Hr'*(R^(-1))*z;
