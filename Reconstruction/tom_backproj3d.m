function tom_backproj3d(volume,image, angle_phi, angle_the, offset, mask)
%TOM_BACKPROJ3D performs a 3D backprojection
%
%   tom_backproj3d(volume,image, angle_phi, angle_the, offset, mask)
%
%PARAMETERS
%
%  INPUT
%   volume              is a 3D volume (single)
%   image               is a 2D image (single)
%   angle_phi           is the projection angle phi
%   angle_the           is the projection angle the
%                       convention see Bronstein
%   offset              is a vector of length three
%  
%  OUTPUT
%
%EXAMPLE
%   for a single-axis reconstruction whith already weighted,
%    aligned projections:
%
%               vol=zeros(1024,1024,512,'single'); 
%               for n=1:52
%                   proj=tom_emread(['TEMP_BPP_' num2str(n) '.em']);
%                   tom_backproj3d(vol,single(proj.Value),0,proj.Header.Tiltangle,[0 0 0]);
%               end;
%
%EXAMPLE
%   for a dual-axis reconstruction whith already weighted,
%    aligned projections:
%
%               vol2=zeros(256,256,128,'single');
%               for n=1:34;
%                       proj=tom_emread(['TEMP_DUAL_' num2str(n) '.em']);
%                       tom_backproj3d(vol2,single(tom_bin(proj.Value,1)),proj.Header.Tiltaxis,proj.Header.Tiltangle,[0 0 0]);
%               end;
%
%
%   tom_backproj3d(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by SN 08/08/02
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

if isa(volume,'single')==0
        error('tom_backproj3d: volume must be single type. No backprojection done!');
end;
if isa(image,'single')==0
        error('tom_backproj3d: projection must be single type. No backprojection done!');
end;
if nargin < 5
     disp('tom_backproj3d: Not enough input parameters. No backprojection done!');
 end
if nargin == 5
     tom_backproj3dc(single(volume),single(image), angle_phi, angle_the, offset);
 end
if nargin == 6
     tom_backproj3dmaskc(single(volume),single(image), angle_phi, angle_the, offset,mask);
end
  
 
