function tom_backproj3d_euler(volume,image,angle_phi,angle_psi,angle_the,offset,mask)
%TOM_BACKPROJ3D_EULER performs a 3D backprojection
%
%   tom_backproj3d_euler(volume,image,angle_phi,angle_psi,angle_the,offset,mask)
%
%   Works for any rotation phi phi theta and z-projection
%
%PARAMETERS
%
%  INPUT
%   volume              is a 3D volume (single)
%   image               is a 2D image (single)
%   angle_phi           is the projection angle phi
%   angle_psi           ...
%   angle_the           is the projection angle the
%                       convention see Bronstein
%   offset              is a vector of length three
%   mask                ...
%  
%  OUTPUT
%
%NOTE 
%   Rot phi psi theta and proj in z-direction is the same as proj psi=90-phi and rot phi + psi
%   inverse: Rot -phi-psi-90+phi ==> -psi-90  ; bakpro= 90-phi psi  
%
%EXAMPLE
%   for a single-axis reconstruction whith already weighted,
%    aligned projections:
%
%               vol=zeros(1024,1024,512,'single'); 
%               for n=1:52
%                   proj=tom_emread(['TEMP_BPP_' num2str(n) '.em']);
%                   tom_backproj3d(vol,proj.Value,phi,psi,theta,[0 0 0]);
%               end;
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

%rotate the projection
image=single(image);
image=single(tom_rotate(image,-angle_psi-90,'linear'));
%image=single(tom_rotate(image,-angle_psi-90,'linear'));



if isa(volume,'single')==0
        error('tom_backproj3d: volume must be single type. No backprojection done!');
end;
if isa(image,'single')==0
        error('tom_backproj3d: projection must be single type. No backprojection done!');
end;
if nargin < 6
     disp('tom_backproj3d: Not enough input parameters. No backprojection done!');
 end
if nargin == 6
     tom_backproj3dc(single(volume),single(image), 90-angle_phi, angle_the, offset);
 end
% if nargin == 7
%      tom_backproj3dmaskc(single(volume),single(image), 90-angle_phi, angle_the, offset,mask);
% end
  