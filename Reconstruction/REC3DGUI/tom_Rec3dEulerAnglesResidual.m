function [EulerAnglesResidual,MaximumResidual,ResidualSpheres,AverageResidualSphere] = ...
tom_Rec3dEulerAnglesResidual(m3d1,m3d2,NumOfMarkers,RotMatrix,ProtocolMode)
%TOM_REC3DEULERANGLESRESIDUAL is a module of TOM_REC3DGUI.
%
%   [EulerAnglesResidual,MaximumResidual,ResidualSpheres,AverageResidualSphere] = ...
%   tom_Rec3dEulerAnglesResidual(m3d1,m3d2,NumOfMarkers,RotMatrix,ProtocolMode)
%
%   It calculates the residual of 3d-rotation and translation 
%   between two markersets
%
%PARAMETERS
%
%  INPUT
%   m3d1                    markerset m3d1
%   m3d2                    markerset m3d2
%   NumOfMarkers            Number of markers
%   RotMatrix               rotation matrix
%   ProtocolMode            'on' or 'off'
%  
%  OUTPUT
%   EulerAnglesResidual     residual matrix of rotation and translation
%   MaximumResidual     	maximum value of residual matrix
%   ResidualSpheres         euclidian length of residuals
%   AverageResidualSphere   mean of euclidian length of residuals
%
%EXAMPLE
%
%REFERENCES
%
%SEE ALSO
%
%   created by ME 01/01/07
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


% Calculate EulerAnglesResidual
for k=1:NumOfMarkers
    
     a = RotMatrix*[m3d1(1,k);m3d1(2,k);m3d1(3,k)];
    
         m3d1rot(1,k) = a(1);
         m3d1rot(2,k) = a(2);
         m3d1rot(3,k) = a(3);

end

for k=1:NumOfMarkers
    
     a = [m3d2(1,k);m3d2(2,k);m3d2(3,k)] - ...
         [m3d1rot(1,k);m3d1rot(2,k);m3d1rot(3,k)];
    
         EulerAnglesResidual(1,k) = a(1);
         EulerAnglesResidual(2,k) = a(2);
         EulerAnglesResidual(3,k) = a(3);

end



% EulerAnglesResidual in absolut values
EulerAnglesResidual = abs(EulerAnglesResidual);

% MaximumResidual
MaximumResidual = max(max(EulerAnglesResidual));

% ResidualSpheres
for k=1:NumOfMarkers
    ResidualSpheres(k) = norm(EulerAnglesResidual(:,k));
end

% AverageResidualSphere
AverageResidualSphere = sum(ResidualSpheres)/NumOfMarkers;



% Display
if (strcmp(ProtocolMode,'on')) == 1
disp('Euler Angles Residual');
disp(EulerAnglesResidual);
disp('Maximum Residual');
disp(MaximumResidual);
disp('Residual Spheres');
disp(ResidualSpheres);
disp('Average Residual Sphere');
disp(AverageResidualSphere);
end

