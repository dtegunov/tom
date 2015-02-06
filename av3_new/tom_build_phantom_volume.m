function vol=tom_build_phantom_volume(vol,tmplate,pos,ang)
%TOM_BUILD_PHANTOM_VOLUME creates ...
%
%   vol=tom_build_phantom_volume(vol,tmplate,pos,ang)
%
%PARAMETERS
%
%  INPUT
%   vol                 ...
%   tmplate             ...
%   pos                 ...
%   ang                 ...
%  
%  OUTPUT
%   vol                 ...
%
%EXAMPLE
%   ... = tom_build_phantom_volume(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by .. mm/dd/yy
%   updated by ..
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


for i=1:size(pos,1)
    act_tmp=tom_rotate(tmplate,[ang(i,:)],'linear');
    mask = ones(size(tmplate)); 
    mask = tom_spheremask(mask, 15, 1, [round(size(act_tmp)./2)]);
    %mask = mask .* 0.5;
    act_tmp = act_tmp .* mask;

    
    vol=tom_paste(vol,act_tmp,pos(i,:));
end;



