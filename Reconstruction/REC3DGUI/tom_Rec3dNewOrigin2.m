function [NewOrigin2] = tom_Rec3dNewOrigin2...
(Origin1,Origin2,Imdim,RotMatrix,ProtocolMode)
%TOM_REC3DNEWORIGIN2 is a module of TOM_REC3DGUI.
%
%   [NewOrigin2] = tom_Rec3dNewOrigin2...
%   (Origin1,Origin2,Imdim,RotMatrix,ProtocolMode)
%
%   It calculates the new origin for tiltseries2.
%
%PARAMETERS
%
%  INPUT
%   Origin1             coordinates of origin1
%   Origin2             coordinates of origin2
%   Imdim               image dimension
%   RotMatrix           rotation matrix
%   ProtocolMode        'on' or 'off'
%  
%  OUTPUT
%   NewOrigin2          new coordinates of origin2
%
%EXAMPLE
%
%REFERENCES
%
%SEE ALSO
%   TOM_REC3DALIGNMENTENGINE
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


% Translate origin from upper left corner to imagecenter
Origin1t = Origin1' - [Imdim/2+1;Imdim/2+1;Imdim/2+1];
Origin2t = Origin2' - [Imdim/2+1;Imdim/2+1;Imdim/2+1];

% Calculate difference vector
Diff = Origin2t - RotMatrix*Origin1t;

% New origin for tiltseries2
NewOrigin2 = Origin2' - Diff;
NewOrigin2 = NewOrigin2';

% % Display NewOrigin2 and new difference with NewOrigin2
% if (strcmp(ProtocolMode,'on')) == 1
%     disp('New Origin of tiltseries2');
%     disp(NewOrigin2);
% end

