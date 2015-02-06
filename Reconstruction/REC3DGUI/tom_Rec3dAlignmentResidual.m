function [AveragePerProjection,AveragePerMarker] = ...
tom_Rec3dAlignmentResidual(ResidualMatrix,NumOfMarkers,NumOfProj)
%TOM_REC3DALIGNMENTRESIDUAL is a module of TOM_REC3DGUI.
%
%   [AveragePerProjection,AveragePerMarker] = ...
%   tom_Rec3dAlignmentResidual(ResidualMatrix,NumOfMarkers,NumOfProj)
%
%   It evaluates alignment residual matrix.
%
%PARAMETERS
%
%  INPUT
%   ResidualMatrix           Residual matrix of alignment
%   NumOfMarkers             Number of markers
%   NumOfProj                Number of projections
%  
%  OUTPUT
%   AveragePerProjection	 Residual average per projection
%   AveragePerMarker	     Residual average per marker
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


% AveragePerProjection
for k=1:NumOfProj
    
sumup = 0;
zaehler = 0;
    
    for j=1:NumOfMarkers
         
         % Exception (-1)
         if ResidualMatrix(j,k) == (-1) 
              zaehler = zaehler + 1;
         else
              sumup = ResidualMatrix(j,k) + sumup;
         end
         
    end
    
    if (NumOfMarkers-zaehler) == 0
         AveragePerProjection(k) = (-1);
    else
         AveragePerProjection(k) = sumup./(NumOfMarkers-zaehler);
    end
    
end


% AveragePerMarker
for k=1:NumOfMarkers

sumup = 0;
zaehler = 0;
    
    for j=1:NumOfProj
        
         % Exception (-1)
         if ResidualMatrix(k,j) == (-1)
              zaehler = zaehler + 1;
         else
              sumup = ResidualMatrix(k,j) + sumup;
         end
    
    end
    
         if (NumOfProj-zaehler) == 0
              AveragePerMarker(k) = (-1);
         else
              AveragePerMarker(k) = sumup./(NumOfProj-zaehler);
         end
         
end

