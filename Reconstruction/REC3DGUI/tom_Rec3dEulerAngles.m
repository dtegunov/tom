function [Psi, Theta, Phi, RotMatrix] = tom_Rec3dEulerAngles...
(m3d1,m3d2,NumOfMarkers,ProtocolMode)
%TOM_REC3DEULERANGLES is a module of TOM_REC3DGUI.
%
%   [Psi,Theta,Phi,RotMatrix] = tom_Rec3dEulerAngles...
%   (m3d1,m3d2,NumOfMarkers,ProtocolMode)
%
%   It calculates Euler angles (Psi,Theta,Phi) and rotation matrix (RotMatrix) 
%   between two markersets of two tiltseries.
%
%PARAMETERS
%
%  INPUT
%   m3d1              markerset m3d1
%   m3d2              markerset m3d2
%   NumOfMarkers      Number of markers
%   ProtocolMode      'on' or 'off'
%  
%  OUTPUT
%   Psi               Euler angle psi
%   Theta             Euler angle theta
%   Phi               Euler angle phi
%   RotMatrix         Rotation matrix
%
%EXAMPLE
%
%REFERENCES
%
%SEE ALSO
%   TOM_REC3DEULERANGLES   TOM_REC3DEULERANGLESRESIDUAL
%   TOM_REC3DNEWORIGIN2
%   TOM_REC3DRIGIDBODYALIGNMENT   TOM_REC3DRIGIDBODYRESIDUAL
%   TOM_REC3DALIGNMENTRESIDUAL
%
%   01/01/07 ME
%
%   Nickell et al., 'TOM software toolbox: acquisition and analysis for electron tomography',
%   Journal of Structural Biology, 149 (2005), 227-234.
%
%   Copyright (c) 2004-2007
%   TOM toolbox for Electron Tomography
%   Max-Planck-Institute of Biochemistry
%   Dept. Molecular Structural Biology
%   82152 Martinsried, Germany
%   http://www.biochem.mpg.de/tom


%--------------------------------------------------------------------------                                              
% CALCULATE QUATERNION MATRIX
%--------------------------------------------------------------------------
A = funA(NumOfMarkers,m3d1,m3d2);
B = funB(NumOfMarkers,m3d1,m3d2);
C = funC(NumOfMarkers,m3d1,m3d2);
D = funD(NumOfMarkers,m3d1,m3d2);
E = funE(NumOfMarkers,m3d1,m3d2);
F = funF(NumOfMarkers,m3d1,m3d2);
G = funG(NumOfMarkers,m3d1,m3d2);
H = funH(NumOfMarkers,m3d1,m3d2);
I = funI(NumOfMarkers,m3d1,m3d2);


Q = [A+B+C   -H+I    F-G    -D+E; ...
     -H+I   A-B-C    D+E     F+G; ...
      F-G    D+E   -A+B-C    H+I; ...
     -D+E    F+G     H+I   -A-B+C];

% Calculate eigenvectors and eigenvalues
[eig_vec_Q,eig_val_Q] = eig(Q);

% Find most positive eigenvalue
eig_val_Q_max = eig_val_Q(1,1);
index_eig_val_Q_max = 1;

for i=2:4
    
    if (eig_val_Q(i,i) >= eig_val_Q_max)
            eig_val_Q_max = eig_val_Q(i,i);
            index_eig_val_Q_max = i;
    end
end

% Find corresponding eigenvector
eig_vec_Q_max = eig_vec_Q(1:4,index_eig_val_Q_max);

q0 = eig_vec_Q_max(1,1);
qx = eig_vec_Q_max(2,1);
qy = eig_vec_Q_max(3,1);
qz = eig_vec_Q_max(4,1);
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------                                              
% CALCULATE ROTATION MATRIX
%--------------------------------------------------------------------------
RotMatrix = [q0*q0+qx*qx-qy*qy-qz*qz  2*(qx*qy+q0*qz)       2*(qx*qz-q0*qy); ...
               2*(qx*qy-q0*qz)     q0*q0-qx*qx+qy*qy-qz*qz  2*(qy*qz+q0*qx); ...
               2*(qx*qz+q0*qy)        2*(qy*qz-q0*qx)  q0*q0-qx*qx-qy*qy+qz*qz];
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------                                              
% CALCULATE EULER ANGLES
%--------------------------------------------------------------------------
Psi   = (180/pi)*atan2(RotMatrix(1,3),-RotMatrix(2,3));
Theta = (180/pi)*acos(RotMatrix(3,3)); % Theta per definition within (0...180) deg
Phi   = (180/pi)*atan2(RotMatrix(3,1),RotMatrix(3,2));

% Special case: R(3,3) is near or is +1, then RotMatrix degenerates 
% to a matrix representing a in-plane rotation around Phi
if -(RotMatrix(3,3)-1)<10e-8
    
    Theta = 0;
    Psi = 0;
    Phi = (180/pi)*atan2(RotMatrix(2,1),RotMatrix(1,1));
    
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Display Euler angles and rotation matrix
%--------------------------------------------------------------------------
if (strcmp(ProtocolMode,'on')) == 1
disp('Euler angles');
disp(['Psi ' num2str(Psi)]);
disp(['Theta ' num2str(Theta)]);
disp(['Phi ' num2str(Phi)]);
disp(' ');
disp('Rot Matrix between tiltseries 1 and tiltseries 2');
disp(RotMatrix);
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------                                              
% SUBFUNCTIONS QUATERNION MATRIX
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
function valA = funA(NumOfMarkers,m3d1,m3d2)

sumA = 0;

for k=1:NumOfMarkers
    
    sumA = sumA + (m3d1(1,k))*(m3d2(1,k));

end

valA = sumA;
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
function valB = funB(NumOfMarkers,m3d1,m3d2)

sumB = 0;

for k=1:NumOfMarkers
    
    sumB = sumB + (m3d1(2,k))*(m3d2(2,k));

end

valB = sumB;
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
function valC = funC(NumOfMarkers,m3d1,m3d2)

sumC = 0;

for k=1:NumOfMarkers
    
    sumC = sumC + (m3d1(3,k))*(m3d2(3,k));

end

valC = sumC;
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
function valD = funD(NumOfMarkers,m3d1,m3d2)

sumD = 0;

for k=1:NumOfMarkers
    
    sumD = sumD + (m3d1(1,k))*(m3d2(2,k));

end

valD = sumD;
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
function valE = funE(NumOfMarkers,m3d1,m3d2)

sumE = 0;

for k=1:NumOfMarkers
    
    sumE = sumE + (m3d1(2,k))*(m3d2(1,k));

end

valE = sumE;
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
function valF = funF(NumOfMarkers,m3d1,m3d2)

sumF = 0;

for k=1:NumOfMarkers
    
    sumF = sumF + (m3d1(1,k))*(m3d2(3,k));

end

valF = sumF;
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
function valG = funG(NumOfMarkers,m3d1,m3d2)

sumG = 0;

for k=1:NumOfMarkers
    
    sumG = sumG + (m3d1(3,k))*(m3d2(1,k));

end

valG = sumG;
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
function valH = funH(NumOfMarkers,m3d1,m3d2)

sumH = 0;

for k=1:NumOfMarkers
    
    sumH = sumH + (m3d1(2,k))*(m3d2(3,k));

end

valH = sumH;
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
function valI = funI(NumOfMarkers,m3d1,m3d2)

sumI = 0;

for k=1:NumOfMarkers
    
    sumI = sumI + (m3d1(3,k))*(m3d2(2,k));

end

valI = sumI;
%--------------------------------------------------------------------------

