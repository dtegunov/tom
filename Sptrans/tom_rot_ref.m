function vol=tom_rot_ref(vol_in,image,im_size_z,angle_range,angle_inkre);
%TOM_ROT_REF creates ...
%
%   vol=tom_rot_ref(vol_in,image,im_size_z,angle_range,angle_inkre)
%
%PARAMETERS
%
%  INPUT
%   vol_in              ...
%   image               ...
%   im_size_z           ...
%   angle_range         ...
%   angle_inkre         ...
%  
%  OUTPUT
%   vol                 ...
%
%EXAMPLE
%   ... = tom_rot_ref(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by ... (author date)
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


middle=floor(size(vol_in,3)./2)+1;
start=middle-round(im_size_z./2);
stop=middle+round(im_size_z./2);
for i=start:stop
    vol_in(:,:,i)=image;
end;

for i_angle=(angle_range(1)+angle_inkre):angle_inkre:angle_range(2)
    i_angle    
    rot_vol=tom_rotate(vol_in,[0 0 i_angle],'linear');
    vol_in=tom_paste(vol_in,rot_vol,[1 1 1],'max');

end;
vol=vol_in;
