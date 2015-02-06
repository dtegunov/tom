function tom_emstack2spiderseries(source_em_file,particle_name,sel_file_name,norm,invert,bin,start)
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
%   norm               gradient/statistics 
%   invert             0 or 1
%   binning            0, 1, ...
%   start              particle starting number
%
%
%EXAMPLE
%   tom_emstack2spiderseries('side.em','parts_','26S.sel')
%   tom_emstack2spiderseries('side.em','parts_','26S.sel','gradient&mean0+1std',1,2)
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
if nargin<4
    norm='no_norm';
end;
if nargin<5
    invert=0;
end;
if nargin<6
    bin=0;
end;

if nargin<7
    start=1;
end;


try
    header=tom_reademheader(source_em_file);
catch
    error('Cannot open file');
end;

if ~tom_isemfile(source_em_file)
    error([source_em_file ' is not an EM-file']);
end;

fp=fopen(sel_file_name,'wt');

disp('ok!');

out_z=start;


for iz=1:header.Header.Size(3)
    
    in=tom_emreadc(source_em_file,'subregion',[1 1 iz],[header.Header.Size(1)-1 header.Header.Size(2)-1 0]);
    in.Value=tom_bin(in.Value,bin);
    
   % in.Value=imresize(in.Value,[64 64]);
    
    if (strcmp(norm,'gradient')||strcmp(norm,'gradient&mean0+1std'))
        in.Value=tom_xmipp_normalize(in.Value,'Ramp');
    end;
    if (invert==1)
        in.Value=-in.Value;
    end;
    if (strcmp(norm,'mean0+1std')||strcmp(norm,'gradient&mean0+1std'))
        in.Value=tom_norm(in.Value,'mean0+1std');
    end;
    
    tom_spiderwrite([particle_name num2str(out_z) '.spi'],in.Value);
    
    fprintf(fp,[particle_name num2str(out_z) '.spi 1 \n']);
    
    out_z=out_z+1;
    
end;
fclose(fp);