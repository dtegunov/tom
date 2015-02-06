function tom_av2_xmipp_comp_frames(part_list,rule,targetFold,targetSel,norm)
%  tom_av2_xmipp_comp_frames combines frames of particles according 2 given
%  rule
%  
%    tom_av2_xmipp_comp_frames(part_list,rule,targetFold)
%  
%  PARAMETERS
%  
%    INPUT
%     part_list                textfile with frame folders
%     rule                     vector describing which frames are combinde
%                              (use 0 to discard)
%     targetFold               folder for combinde frames
%     targetSel                sel-file for the new frames
%     norm                     (mean0+1std-bg)  norming 4 the sumed frames  
%  
%  EXAMPLE
%  
%    unix('cat sum.sel | awk ''{gsub(".spi","-spi"); print $1}''  > list.txt');
%    
%    zz=1; for i=1:13; for ii=1:3; rule(zz)=i;zz=zz+1; end; end;
%    rule(zz)=0; %cluster 3 and discard the last !
%    tom_av2_xmipp_comp_frames('list.txt',rule,'pFramesCl3','pcl3.sel');
%
%  REFERENCES
%  
%  SEE ALSO
%     ...
%  
%     created by FB 25/04/13
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

if (nargin < 5)
    norm='mean0+1std-bg';
end;



warning off; mkdir(targetFold); warning on;

list=importdata(part_list);
[aa bb cc]=fileparts(list{1});
p_tmp_name=strtok(bb,'-');
im_tmp=tom_spiderread([list{1} filesep p_tmp_name '_F_1.spi']);
sz_im=size(im_tmp.Value); 

if (findstr(norm,'bg'))
    mask=tom_spheremask(ones(sz_im),round(sz_im(1)./2)-1)==0;
    norm=strrep(norm,'-bg','');
else
    mask=ones(sz_im);
end;

u_rule=unique(rule); 
idx=ismember(u_rule,0);
u_rule=u_rule(find(idx==0));


disp('Clustering Frames to rule:');
disp(rule);

disp(['Wrinting 2 : ' targetFold]);
disp(['Wrinting 2 : ' targetSel]);
disp(' ');
fid=fopen(targetSel,'wt');

for i=1:length(list)
    [a b c]=fileparts(list{i});
    new_part_fold=[targetFold filesep b];
    warning off; mkdir(new_part_fold); warning on;
    [aa bb cc]=fileparts(new_part_fold);
    p_tmp_name=strtok(bb,'-');
    fra_base=[new_part_fold filesep p_tmp_name '_F_'];
    disp(['Processing ' list{i}]);    
    for ii=1:length(u_rule)
        fr_idx=find(u_rule(ii)==rule);
        buffer=zeros(sz_im(1),sz_im(2),length(fr_idx));
        for iii=1:length(fr_idx)
            [aa bb cc]=fileparts(list{i});
            p_tmp_name=strtok(bb,'-');
            frame_name=[list{i} filesep p_tmp_name '_F_' num2str(fr_idx(iii)) '.spi'];
            tmp=tom_spiderread(frame_name);
            buffer(:,:,iii)=tmp.Value;
        end;
        s_tmp=sum(buffer,3)./size(buffer,3);
        fr_cl=tom_norm(s_tmp,norm,mask);
        tom_spiderwrite([fra_base num2str(ii) '.spi'],fr_cl);
        fprintf(fid,'%s 1\n',[fra_base num2str(ii) '.spi']);
    end;
end;

fclose(fid);

















