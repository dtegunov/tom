function [projection]=tom_proj3d(volume,angle)
%TOM_PROJ3D calculates 2D projection(s) of a 3D volume
%
%   [projection]=tom_proj3d(volume,angle)
%
%PARAMETERS
%
%  INPUT
%   volume              3D volume
%   angle               [phi the] projection direction
%                        angles are in degree
%  
%  OUTPUT
%   projection          projection of the volume is calculated after
%                        rotation
%
%    Note that the command
%        tom_proj3d at angle = [phi the]
%    is equivalent to
%        tom_rotate at angle = [90-phi phi-90 the]
%        tom_proj3d at angle = [0 0]
%
% the 'inverse' operation to the backprojection operator such as:
%      tom_backproj3d(volume,proj,proj_tiltaxis,proj_tiltangle,[0 0 0]);
% is:
%      volume_rot=tom_rotate(volume,[-proj_tiltaxis 0 0]);
%      new_proj=tom_proj3d(volume_rot,[proj_tiltaxis proj_tiltangle]);
%
%EXAMPLE
%   [projection]=tom_proj3d(ones(64,64,64,'single'),[20 30])
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by SN 04/08/05
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
        error('tom_proj3d: volume must be single type. No backprojection done!');
end;
for i=1:size(angle,1)
    projection(:,:,i)=sum(tom_rotate(volume,[90-angle(i,1) angle(i,1)-90 angle(i,2)]),3);
end;
