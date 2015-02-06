function tom_fit_ctf_folder2list(basepath,output_textfile)
%TOM_FIT_CTF_FOLDER2LISt converts *.mat files 2 ascii list
%
%
%   tom_fit_ctf_folder2list(basepath,output_textfile)
%
%PARAMETERS
%
%  INPUT
%   basepath             basepath of the mat-files containing teh
%   output_textfile      output textfile
% 
%
%EXAMPLE
%   
%  tom_fit_ctf_folder2list('/fs/sandy02/lv08/pool/pool-tpp/TPP_tom_acquisition/TPP_240709_clean_2/792009/low/*.mat','../out.txt');
%
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

fp=fopen(output_textfile,'w');

[a b c]=fileparts(basepath);

pp=dir(basepath);


for i=1:size(pp,1)
    load([a '/' pp(i).name]);
    try
        fprintf(fp,[a '/' strrep(pp(i).name,'.mat','') ' ' num2str(st_out.sel.selected) '\n']);
    catch ME
        disp(ME.message);
    end;
    
end;

fclose(fp);