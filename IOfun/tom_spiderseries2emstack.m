function tom_spiderseries2emstack(sel_file_name,output_em_file)
%TOM_EMSTACK2SPIDERSERIES converts a Spider series defined in a .sel file
%into an EM stack
%
%   tom_spiderseries2emstack(sel_file_name,output_em_file)
%
%PARAMETERS
%
%  INPUT
%   sel_file_name     filename of selfile 
%  OUTPUT
%   output_em_file    filename of new EM stack
%
%EXAMPLE
%    tom_spiderseries2emstack('./data/26S.sel','w.em')
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

i=0;
name='o';
try
    fp=fopen(sel_file_name,'rt');
catch
    error('Cannot open sel-file');
end;
while ~isequal(name,'')
try
    [name c]= fscanf(fp,'%s\n',1);
    [nr c]= fscanf(fp,'%s\n',1);
    i=i+1;
catch
    fclose(fp);
    break;
end;
end;
fp=fopen(sel_file_name,'rt');
[name c]= fscanf(fp,'%s\n',1);
in=tom_spiderread(name);
fclose(fp);


tom_emwritec(output_em_file,[size(in.Value,1) size(in.Value,2) i-1],'new','single');
fp=fopen(sel_file_name,'rt');
for idx=1:i-1
    if idx./1000==round(idx./1000) disp([num2str(idx) ' particles done']);end;
    [name c]= fscanf(fp,'%s\n',1);
    [nr c]= fscanf(fp,'%s\n',1);
    in=tom_spiderread(name);
    tom_emwritec(output_em_file,in.Value,'subregion',[1 1 idx],[size(in.Value,1) size(in.Value,2) 1]);
    
end;
fclose(fp);