function [ResidualMatrix, WarpAlignment] = tom_Rec3dRigidBodyResidual...
(MarkersOnProj,NumOfProj,NumOfMarkers,Origin,Imdim,Tiltaxis,Tiltangles,m3d,tx,ty)
%TOM_REC3DRIGIDBODYRESIDUAL is a module of TOM_REC3DGUI.
%
%   [ResidualMatrix,WarpAlignment] = tom_Rec3dRigidBodyResidual...
%   (MarkersOnProj,NumOfProj,NumOfMarkers,Origin,Imdim,Tiltaxis,Tiltangles,m3d,tx,ty)
%
%   It calculates residual matrix of a rigid body alignment.
%
%PARAMETERS
%
%  INPUT
%   MarkersOnProj     xy-coordinates of markers measured in tom_setmark of each tiltseries
%   NumOfProj         Number of projections of each tiltseries
%   NumOfMarkers      Number of markers
%   Origin            assigned coordinate origin of each tiltseries
%   Imdim             dimension of projections
%   Tiltaxis          calculated angle between projection and specimen frame around beam axis
%   Tiltangles        projection tilting angle 
%   m3d               calculated 3d-coordinates of markers in specimen frame
%   tx                calculated shift in x-direction
%   ty                calculated shift in x-direction
%  
%  OUTPUT
%   ResidualMatrix	  residuals per projection and per marker
%   WarpAlignment	  xy-coordinates of markers calculated with alignment
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


% Angles in radiant
Tiltaxis = (pi./180).*Tiltaxis;
Tiltangles = (pi./180).*Tiltangles;


% Calculate WarpAlignment
for j=1:NumOfMarkers
    for i=1:NumOfProj
         
         % Exception (-1)
         if MarkersOnProj(j,i,1) == (-1) || MarkersOnProj(j,i,2) == (-1)
              
              WarpAlignment(j,i,1) = (-1);
              WarpAlignment(j,i,2) = (-1);
         
         else
             
         % Rigid body rotation and projection operator
         P = [(sin(Tiltaxis(i))).^2.*cos(Tiltangles(i))+(cos(Tiltaxis(i))).^2 ...
                  sin(Tiltaxis(i)).*cos(Tiltaxis(i)).*(1-cos(Tiltangles(i))) ...
                  sin(Tiltaxis(i)).*sin(Tiltangles(i)); ...
                  sin(Tiltaxis(i)).*cos(Tiltaxis(i)).*(1-cos(Tiltangles(i))) ...
                 (cos(Tiltaxis(i))).^2.*cos(Tiltangles(i))+(sin(Tiltaxis(i))).^2 ...
                -cos(Tiltaxis(i)).*sin(Tiltangles(i))];
         
         % Calculated marker coordinates
         a = P*[m3d(1,j)+Origin(1)-(Imdim./2+1);...
                    m3d(2,j)+Origin(2)-(Imdim./2+1);...
                    m3d(3,j)+Origin(3)-(Imdim./2+1)] + ...
                   [tx(i)+(Imdim./2+1);ty(i)+(Imdim./2+1)];
         
         WarpAlignment(j,i,1) = a(1);
         WarpAlignment(j,i,2) = a(2);
         
         end
         
    end
end


% Calculate ResidualMatrix
for j=1:NumOfMarkers
    for i=1:NumOfProj
    
         % Exception (-1)
         if MarkersOnProj(j,i,1) == (-1) || MarkersOnProj(j,i,2) == (-1)
              
              ResidualMatrix(j,i) = (-1);
         
         else
              
              ResidualMatrix(j,i) = sqrt((MarkersOnProj(j,i,1) - WarpAlignment(j,i,1)).^2 + ...
                                                        (MarkersOnProj(j,i,2) - WarpAlignment(j,i,2)).^2);
         end

    end
end

