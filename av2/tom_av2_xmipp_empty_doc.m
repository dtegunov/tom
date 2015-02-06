function st=tom_av2_xmipp_empty_doc(sel_file_name,output_doc)
%  tom_av2_xmipp_empty_doc creates an empty doc file ...speed up
%  xmipp_projmatch header extract can be skippt 
%    
%       tom_av2_xmipp_empty_doc(sel_file_name,output_doc)
%    
%    PARAMETERS
%    
%      INPUT
%       sel_file_name   sel file name
%       output_doc      output doc name          
%      
%      OUTPUT
%       st             empty doc struct 
%
%    
%    EXAMPLE
%        tom_av2_xmipp_empty_doc('my.sel','out_doc.doc');
%  
%        %gen empty doc 4 programming
%        tmp_doc=tom_av2_xmipp_empty_doc(100);
%
%   REFERENCES
%    
%    NOTE:
%   
%    
%    
%  
%    SEE ALSO
%       tom_xmippdocread,tom_spiderread
%    
%       created by FB 15/02/10
%    
%       Nickell et al., 'TOM software toolbox: acquisition and analysis for electron tomography',
%       Journal of Structural Biology, 149 (2005), 227-234.
%    
%       Copyright (c) 2004-2007
%       TOM toolbox for Electron Tomography
%       Max-Planck-Institute of Biochemistry
%       Dept. Molecular Structural Biology
%       82152 Martinsried, Germany
%       http://www.biochem.mpg.de/tom

if (nargin < 2)
    output_doc='';
end;


if (isnumeric(sel_file_name))
    for i=1:sel_file_name
        sel.textdata{i}='';
    end;
else
    %read sel
    sel=importdata(sel_file_name);
end;

%head_st=get_header_info('maxCC');
st=alloc_st(length(sel.textdata),{'maxCC'});

head=' ; Headerinfo columns: rot (1), tilt (2), psi (3), Xoff (4), Yoff (5), Ref (6), Flip (7), maxCC (8)';

for i=1:length(st);
    st(i).name=sel.textdata{i};
end;
st(1).header=head;

if (isempty(output_doc)==0)
    tom_xmippdocwrite(output_doc,st);
end;



function st=alloc_st(num_of_entry,header_st)

flip=cell(num_of_entry,1);
ref=cell(num_of_entry,1);
yoff=cell(num_of_entry,1);
xoff=cell(num_of_entry,1);
psi=cell(num_of_entry,1);
tilt=cell(num_of_entry,1);
rot=cell(num_of_entry,1);
part_idx=cell(num_of_entry,1);
name=cell(num_of_entry,1);
run_num=cell(num_of_entry,1);
cols=cell(num_of_entry,1);
tmp_1=cell(num_of_entry,1);

for i=1:length(header_st)
    header_st{i}=strrep(header_st{i},'/','_');
end;


for i=1:num_of_entry
    flip{i}=0;
    ref{i}=0;
    yoff{i}=0;
    xoff{i}=0;
    psi{i}=0;
    tilt{i}=0;
    rot{i}=0;
    part_idx{i}=0;
    tmp_1{i}=0;
    name{i}='                                                                                                                                                                                                             ';
    cols{i}=8;
    run_num{i}=i;
end;

if (length(header_st)==1)
    st=struct('flip',flip,'ref',ref,'yoff',yoff,'xoff',xoff,'psi',psi,'tilt',tilt,'rot',rot,'part_idx',part_idx,'name',name,'cols',cols,'run_num',run_num,header_st{1},tmp_1);
end;

if (length(header_st)==2)
    st=struct('flip',flip,'ref',ref,'yoff',yoff,'xoff',xoff,'psi',psi,'tilt',tilt,'rot',rot,'part_idx',part_idx,'name',name,'cols',cols,'run_num',run_num,header_st{1},tmp_1,header_st{2},tmp_1);
end;

if (length(header_st)==3)
    st=struct('flip',flip,'ref',ref,'yoff',yoff,'xoff',xoff,'psi',psi,'tilt',tilt,'rot',rot,'part_idx',part_idx,'name',name,'cols',cols,'run_num',run_num,header_st{1},tmp_1,header_st{2},tmp_1,header_st{3},tmp_1);
end;
 
if (length(header_st)==4)
    st=struct('flip',flip,'ref',ref,'yoff',yoff,'xoff',xoff,'psi',psi,'tilt',tilt,'rot',rot,'part_idx',part_idx,'name',name,'cols',cols,'run_num',run_num,header_st{1},tmp_1,header_st{2},tmp_1,header_st{3},tmp_1,header_st{4},tmp_1);
end

if (length(header_st)==5)
    st=struct('flip',flip,'ref',ref,'yoff',yoff,'xoff',xoff,'psi',psi,'tilt',tilt,'rot',rot,'part_idx',part_idx,'name',name,'cols',cols,'run_num',run_num,header_st{1},tmp_1,header_st{2},tmp_1,header_st{3},tmp_1,header_st{4},tmp_1,header_st{5},tmp_1);
end






function head_st=get_header_info(line_in)


idx=strfind(line_in,',');

if (strfind(line_in(idx(6):idx(7)),'Flip')==0)
    disp('strange Header in doc !!');
    return;
end;
zz=0;
for i=7:length(idx)
    if (i==length(idx))
        tmp=line_in(idx(i):end);
    else
        tmp=line_in(idx(i):idx(i+1));
    end;
    idx_tmp=strfind(tmp,' ');
    zz=zz+1;
    head_st{zz}=tmp(idx_tmp(1)+1:idx_tmp(2)-1);
end;



