function tom_emstack2emdir(source_em_file,particle_name, sel_file_name)

%TOM_EMSTACK2SPIDERSERIES converts an EM stack to a Spider series of images
%and creates a .sel file
%
%   tom_emstack2spiderseries(source_em_file)
%
%PARAMETERS
%
%  INPUT
%   source_em_file     filename of volume (EM format!)
%   particle_name      filename of individual images (Spider format),
%                       without extension!
%   sel_file_name      filename of sel file with extension!
%  OUTPUT
%   series of images 
%   .sel file
%
%EXAMPLE
%   tom_emstack2spiderseries('side.em','parts_','26S.sel')
%
%REFERENCES
%
%SEE ALSO
%   tom_mrc2emstack, tom_emread, tom_emwrite
%
%   created by SN 01/08/07
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
try
    header=tom_reademheader(source_em_file);
catch
    error('Cannot open file');
end;
fp=fopen(sel_file_name,'wt');
m = tom_sphere([80 80],38);
for iz=1:header.Header.Size(3)

    in=tom_emreadc(source_em_file,'subregion',[1 1 iz],[header.Header.Size(1)-1 header.Header.Size(2)-1 0]);
    
    %in.Value = -single(in.Value);
    %in.Value = tom_xmipp_scale_pyramid(in.Value,'reduce',1);
    %in.Value = tom_xmipp_normalize(in.Value,'ramp',m);
    tom_emwrite([particle_name num2str(iz) '.em'],in.Value);
%    m = mean2(in.Value);
%    s = std2(in.Value);
    
    %tom_emwrite([particle_name num2str(iz) '.em'],tom_norm(tom_limit(in.Value,m-3*s,m+3*s),'mean0+1std'));
   
    fprintf(fp,[particle_name num2str(iz) '.em \n']);
    
end;
fclose(fp);