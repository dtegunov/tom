function st=tom_xmippdocread(filename)
% TOM_XMIPPDOCREAD reads xmipp doc file
%  
%     st=tom_xmippdocread(filename)
%  
%  PARAMETERS
%  
%    INPUT
%     filename      filename of xmipp doc file
%    
%    OUTPUT
%     st           matlab struct 
%                          
%  
%  EXAMPLE
%      st=tom_xmippdocread('model_it000001.doc');
%
%  REFERENCES
%  
%  NOTE:
% 
%  
%  
%
%  SEE ALSO
%     tom_spiderread
%  
%     created by FB 15/02/10
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
% 



try
    use_import=0;
    
    [a b c]=fileparts(filename);
    
    %chop doc in parts
    tmp_out_names=[a b 'x34rt_names'  c 'tmp' ];
    unix(['grep " ; " '  filename ' > ' tmp_out_names]);
    
    if (use_import==1)
        all_names=importdata(tmp_out_names);
    else
        fid=fopen(tmp_out_names);
        for i=1:(281474976710655-1)
            tmp=fgetl(fid);
            if ~ischar(tmp)
                break;
            end;
            all_names{i}=tmp(2:end);
        end;
        fclose(fid);
    end;
    unix(['rm ' tmp_out_names]);
    
    tmp_out_data=[a b 'x34rt_data'  c 'tmp' ];
    unix(['grep -v " ; " '  filename ' > ' tmp_out_data]);
    
    if (use_import==1)
        all_data=importdata(tmp_out_data);
    else
        fid=fopen(tmp_out_data);
        tt=textscan(fid,'%f %f %f %f %f %f %f %f %f %f');
        all_data=[tt{1} tt{2} tt{3} tt{4} tt{5} tt{6} tt{7} tt{8} tt{9} tt{10}];
        fclose(fid);
        unix(['rm ' tmp_out_data]);
    end;
    head_st=get_header_info(all_names{1});
    
    st=alloc_st(size(all_data,1),head_st);
    for i=1:size(all_data,1)
        line=strrep(all_names{i+1},';','');
        line=strtrim(line);
        st(i).run_num=all_data(i,1);
        st(i).cols=all_data(i,2);
        st(i).name=line;
        st(i).flip=all_data(i,9);
        if (isempty(head_st))
            st(i).weight=all_data(i,8);
        else
            st(i).ref=all_data(i,8);
        end;
        st(i).yoff=all_data(i,7);
        st(i).xoff=all_data(i,6);
        st(i).psi=all_data(i,5);
        st(i).tilt=all_data(i,4);
        st(i).rot=all_data(i,3);
        for ii=1:length(head_st)
            st(i)=setfield(st(i),strrep(head_st{ii},'/','_'),all_data(i,ii+9));
        end;
        try
            [a b c]=fileparts(line);
            num=strfind(b,'_');
            part_idx=str2double(b(max(num)+1:end));
            st(i).part_idx=part_idx;
        catch ME
            disp(ME.message);
            disp('Warning Cannot extract part nr !!');
            disp('Use name_1.ext conv!');
        end;
    end;
    st(1).header=all_names{1};
   
catch ME
    
    disp(ME.message);
    if (isempty(a))
        a='./';
    end
    disp(['Cannot write tmp -files in: ' a]);
    disp('Switching 2 mat fuction ...takes longer!!!');
    disp(['use chmod ugo+rwx ' a]);
    [a b]=unix(['wc -l ' filename]);
    b=str2double(strtok(b,' '));
    num_of_entry=(b-1)./2;
    
 
    try
        fid=fopen(filename,'rt');
    catch Me
        error('Cannot open doc-file');
    end;
    
    
    line=fgetl(fid);
    head_st=get_header_info(line);
    st=alloc_st(num_of_entry,head_st);
    st(1).header=line;
    num=1;
    zz=1;
    while 1
        
        line=fgetl(fid);
        
        if line==-1
            break;i
        end;
        
        if  mod(num,2)
            line=strrep(line,';','');
            line=strtrim(line);
            st(zz).name=line;
            try
                [a b c]=fileparts(line);
                nums=strfind(b,'_');
                part_idx=str2double(b(max(nums)+1:end));
                st(zz).part_idx=part_idx;
            catch ME
                disp(ME.message);
                disp('Warning Cannot extract part nr !!');
                disp('Use name_1.ext conv!');
            end;
        else
            vals=sscanf(line,'%f');
            st(zz).run_num=vals(1);
            st(zz).cols=vals(2);
            st(zz).flip=vals(9);
            if (isempty(head_st))
                st(zz).weight=vals(8);
            else
                st(zz).ref=vals(8);
            end;
            st(zz).yoff=vals(7);
            st(zz).xoff=vals(6);
            st(zz).psi=vals(5);
            st(zz).tilt=vals(4);
            st(zz).rot=vals(3);
            for ii=1:length(head_st)
                st(zz)=setfield(st(zz),strrep(head_st{ii},'/','_'),vals(ii+9));
            end;
            zz=zz+1;
        end;
        num=num+1;
    end;
    fclose(fid);
   
end;


 %add header and check for unique part_names!!
    
    if (length([st(:).part_idx])==length(unique([st(:).part_idx])) )
        st(1).part_idx_unique=1;
    else
        disp('Warning part idx is not unique !!');
        st(1).part_idx_unique=0;
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

if (isempty(header_st)==0)
    for i=1:length(header_st)
        header_st{i}=strrep(header_st{i},'/','_');
    end;
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
    cols{i}=0;
    run_num{i}=0;
end;

if (length(header_st)==0)
    st=struct('flip',flip,'weight',ref,'yoff',yoff,'xoff',xoff,'psi',psi,'tilt',tilt,'rot',rot,'part_idx',part_idx,'name',name,'cols',cols,'run_num',run_num);
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

if (length(idx) > 6 )
    if (strfind(line_in(idx(6):idx(7)),'Flip')==0)
        disp('strange Header in doc !!');
        return;
    end;
else
    if (strfind(line_in(idx(5):idx(6)),'Flip')==0)
        disp('strange Header in doc !!');
        return;
    end;
    head_st='';
end;

if length(idx) > 6
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
end;















