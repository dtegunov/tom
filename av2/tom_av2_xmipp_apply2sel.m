function tom_av2_xmipp_apply2sel(sel_file,function_name,find_what,replace_with,output_sel)
%TOM_AV2_XMIPP_APPLY2SEL applies custom function 2 all particles in a cell
%
%   tom_av2_xmipp_apply2sel(sel_file,function_name,find_what,replace_with)
%
%  TOM_AV2_XMIPP_APPLY2SEL applies custom function 2 all particles in a sel
%                          useful for applying a mask 2 all particles in a sel 
%   
%  
%
%PARAMETERS
%
%  INPUT
%   f_sel              filename of the input sel
%   function_name      name of the custom function that should be applied 2 each particle  
%   find_what          string in the sel which should be replaced
%   replace_with       new string (according 2 find and replace in a text editor)
%   output_sel         name of the new sel-file
%       
%
%EXAMPLE
%     
%  
%   
%   
%   
%
%REFERENCES
%
%SEE ALSO
%   tom_aling2d
%
%   created by fb ...ole !!
%
%   Nickell et al., 'TOM software toolbox: acquisition and analysis for electron tomography',
%   Journal of Structural Biology, 149 (2005), 227-234.
%
%   Copyright (c) 2004-2007
%   TOM toolbox for Electron Tomography
%   Max-Planck-Institute of Biochemistry
%   Dept. Molecular Structural Biology
%   82152 Martinsried, Germany
d=importdata(sel_file);
tmp=tom_spiderread(d.textdata{1});  



fid=fopen(new_sel,'wt');
tic;







