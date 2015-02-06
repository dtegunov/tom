function tom_spiderdir2emdir(source_dir,dest_dir)

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

sel_file_name = [dest_dir '/files.txt'];

fp=fopen(sel_file_name,'wt');
for i=1:11286
    in=tom_spiderread([source_dir '/26S_' num2str(i) '.spi']);
    
    tom_emwrite([dest_dir '/26S_' num2str(i) '.em'],in.Value);
    
    fprintf(fp,['26S_' num2str(i) '.em\n']);
    
end;
fclose(fp);