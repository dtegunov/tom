function tom_av2_htl_conv2index(input_htl)
%tom_av2_htl_conv2index converts an old htl format 
%   
%  tom_av2_high2low_transform2(htl_filename)
%
%  TOM_AV2_HTL_CONV2INDEX converts an old htl format 2 colums with the part
%  idx are added (function wroks inplace)
%  
%
%PARAMETERS
%
%  INPUT
%   htl_filename        filename of htl file
%
%  OUTPUT
%
%EXAMPLE
%     
%  
%
%
%REFERENCES
%
%SEE ALSO
%    
%  tom_av2_high2low_match2.m
%
%   created by  fb
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


disp(['reading ' input_htl]);
in_htl=importdata(input_htl);


for i=1:length(in_htl)
    % create high list
    [start rest]=strtok(in_htl{i});
    
    names_h{i}=deblank(start);
    [a b c]=fileparts(names_h{i});
    [a d]=strtok(b,'_');
    names_h_tmp(i)=str2double(strrep(d,'_',''));
    
    % create low list
    tmp_l=deblank(strrep(rest,',',''));
    names_l{i}=tmp_l(2:end);
    [a b c]=fileparts(names_l{i});
    [a d]=strtok(b,'_');
    names_l_tmp(i)=str2double(strrep(d,'_',''));

end;

fp=fopen(input_htl,'wt');
for i=1:length(in_htl)
    fprintf(fp,'%s %s %d %d\n',names_h{i},names_l{i},names_h_tmp(i),names_l_tmp(i));
end;
fclose(fp);






