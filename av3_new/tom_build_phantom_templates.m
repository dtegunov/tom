function tom_build_phantom_templates(tmplate,pos,ang,filename,partsize)
%TOM_BUILD_PHANTOM_TEMPLATES generates rotated and shifted volumes for testing purposes
%
%   tom_build_phantom_templates(tmplate,pos,ang,filename,partsize)
%
%PARAMETERS
%
%  INPUT
%   tmplate             ...
%   pos                 ...
%   ang                 ...
%   filename            ...
%   partsize            ...
%  
%  OUTPUT
%
%EXAMPLE
%   tom_build_phantom_templates(...);
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


tmplate = tom_rescale3d(tmplate,partsize);
for i=1:size(pos,1)

    act_tmp = tom_move(tmplate, pos(i,:));
    act_tmp=tom_rotate(act_tmp,ang(i,:),'linear');
    
%    mask = ones(size(tmplate)); 
%    mask = tom_spheremask(mask, 45, 1, [round(size(act_tmp)./2)]);
    %mask = mask .* 0.5;
    %act_tmp = act_tmp .* mask;
   
    

    tom_emwrite(strcat(filename,num2str(i),'.em'),act_tmp);
    
end;



