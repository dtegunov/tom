function tom_av2_xmipp_move_idx_sel(in_sel,offset,ouputbase)
% TOM_AV2_XMIPP_MOVE_IDX_SEL gives an offset to a sel
%  
%     tom_av2_xmipp_move_idx_sel(in_sel,offset,ouputfold)
%  
%    TOM_AV2_XMIPP_MOVE_IDX_SEL gives an offset to a sel and copies the
%    particles
%    
%    
%  
%  PARAMETERS
%  
%    INPUT
%     in_sel       filename of the input sel
%     offset       string in the cell which should be replaced
%     ouputfold    new string (according 2 find and replace in a text editor)
%                          
%  
%  EXAMPLE
%       
%    
%     
%     
%     
%  
%  REFERENCES
%  
%  SEE ALSO
%     
%  
%     created by fb ...ole !!
%  
%     Nickell et al., 'TOM software toolbox: acquisition and analysis for electron tomography',
%     Journal of Structural Biology, 149 (2005), 227-234.
%  
%     Copyright (c) 2004-2007
%     TOM toolbox for Electron Tomography
%     Max-Planck-Institute of Biochemistry
%     Dept. Molecular Structural Biology
%     82152 Martinsried, Germany
%     http://www.biochem.mpg.de/tom

sel=importdata(in_sel);

fp=fopen([in_sel '.mov'],'wt');
fpl=fopen('log.txt','wt');

for i=1:length(sel.textdata)
    tmp_str=sel.textdata{i};
    [a b c]=fileparts(tmp_str);
    pos=strfind(b,'_');
    index=str2double(b(pos+1:end));
    old_name=tmp_str;
    new_name=[ouputbase num2str(index+offset) '.spi']; 
    unix(['cp ' old_name ' ' new_name]);
    fprintf(fp,'%s 1\n',new_name);
    fprintf(fpl,'%s  %s\n',old_name,new_name);
end;
    
fclose(fp);
fclose(fpl);
