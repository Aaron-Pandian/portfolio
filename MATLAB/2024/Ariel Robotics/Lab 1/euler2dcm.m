function [R_BW] = euler2dcm(e)
% euler2dcm : Converts Euler angles phi = e(1), theta = e(2), and psi = e(3)
%             (in radians) into a direction cosine matrix for a 3-1-2 rotation.
%
% Let the world (W) and body (B) reference frames be initially aligned.  In a
% 3-1-2 order, rotate B away from W by angles psi (yaw, about the body Z
% axis), phi (roll, about the body X axis), and theta (pitch, about the body Y
% axis).  R_BW can then be used to cast a vector expressed in W coordinates as
% a vector in B coordinates: vB = R_BW * vW
%
% INPUTS
%
% e ---------- 3-by-1 vector containing the Euler angles in radians: phi =
%              e(1), theta = e(2), and psi = e(3)
%
%
% OUTPUTS
%
% R_BW ------- 3-by-3 direction cosine matrix 
% 
%+------------------------------------------------------------------------------+
% References: Attitude Transformations. VectorNav. (n.d.). 
% https://www.vectornav.com/resources/inertial-navigation-primer/math-fundamentals/math-attitudetran 
%
% Author: Aaron Pandian
%+==============================================================================+  

phi = e(1,1);
theta = e(2,1);
psi = e(3,1);

% Method 1 
R1e = [1 0 0; 0 cos(phi) sin(phi); 0 -sin(phi) cos(phi)];
R2e = [cos(theta) 0 -sin(theta); 0 1 0; sin(theta) 0 cos(theta)];
R3e = [cos(psi) sin(psi) 0; -sin(psi) cos(psi) 0; 0 0 1];

% Method 2 
% function [R] = rotationMatrix(aHat,phi)
% 
% function [uCross] = crossProductEquivalent(u)
% u1 = u(1,1); % extracted values from row 1, column 1
% u2 = u(2,1);
% u3 = u(3,1);
% uCross = [0 -u3 u2; u3 0 -u1; -u2 u1 0];
% end
% 
% I = [1 0 0; 0 1 0; 0 0 1];
% aHatTranspose = transpose(aHat);
% R1 = cos(phi)*I;
% R2 = (1-cos(phi))*aHat*aHatTranspose;
% % or I + crossProductEquivilant(aHat)^2 = a*a^T
% R3 = sin(phi)*crossProductEquivalent(aHat);
% R = R1+R2-R3;
% end
% 
% R3e = rotationMatrix([0;0;1],psi);
% R1e = rotationMatrix([1;0;0],phi);
% R2e = rotationMatrix([0;1;0],theta);

% Using 3-1-2 Rotation 
R_BW = R2e*R1e*R3e;

% Method 3 with Toolbox
% R_BW = eul2rotm([phi,theta,psi],'XYZ');

end









