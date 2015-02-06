function tom_phantom_proj(infile, fileextension, outfile, filerange, anglerange, snratio)
%TOM_PHANTOM_PROJ creates ...
%
%   tom_phantom_proj(infile, fileextension, outfile, filerange, anglerange, snratio)
%
%PARAMETERS
%
%  INPUT
%   infile              ...
%   fileextension       ...
%   outfile             ...
%   filerange           ...
%   anglerange          ...
%   snratio             ...
%  
%  OUTPUT
%
%EXAMPLE
%   tom_phantom_proj(...);
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


if strcmp(fileextension(1), '.') ~= 1
    fileextension = strcat('.', fileextension);
end;

c = 1;

filter = 1;

filter=round(filter*(28./2));
[w_func,c_mask]=create_masks('',filter,0);

for i=filerange
   
   invol = tom_emreadc(strcat(infile,num2str(i),fileextension));
   outvol = zeros(size(invol.Value),'single');
   
   for j=anglerange
       projection = tom_proj3d(invol.Value, [j 0]);
       %noise = rand(size(projection));
       %projection = projection + (noise .* snratio);
       
       projection=tom_weight3d('analytical',projection,filter,'w_func',w_func,'c_mask',c_mask);  
       projection = single(projection); 
       tom_backproj3d(outvol,projection,  0, j, [0 0 0]); 
   end
   tom_emwrite(strcat(outfile, num2str(i), fileextension), outvol);
   disp([num2str(c) ' of ' num2str(size(filerange,2)) ' done.']);
   c = c + 1;
   
end




function [w_func,c_mask]= create_masks(name,filter,pre_binning)
%pic_h=tom_emreadc(name);
%dim_h=size(pic_h.Value);
dim_h = [28 28];
[xh,yh] = ndgrid(-dim_h(1)./2:((dim_h(2)./2)-1));
w_func = tom_norm(abs(xh),1);
c_mask = sqrt((xh).^2 + (yh).^2) <= filter;