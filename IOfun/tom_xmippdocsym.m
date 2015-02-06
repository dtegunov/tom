function tom_xmippdocsym(input_doc,output_doc,sym,flag,sym_parts_path,cp_flag,htl_filename,dir,sym_parts_path_htl)
%tom_xmippdocsym syms a doc file
%   file
%
%  tom_xmippdocsym(input_doc,output_doc,sym)
%
%  TOM_XMIPPDOCSYM syms a doc file
%                  
%                  
%PARAMETERS
%
%  INPUT
%   input_doc            *.doc filename 
%   output_doc            name of output em stack
%   sym                   symmetry c1..cn
%   flag                  for upright positon ('upright')
%   sym_parts_path        (opt.) path for the sym particels
%   cp_flag               (lns) flag for part duplication cp for copy lns for soft link 
%   htl_filename          (opt.) high2low filename (will be extended 2 new particles) 
%   dir                   (high2low) direction   
%   sym_parts_path_htl    path for particles
%
%  OUTPUT
%
%EXAMPLE
%     
%  tom_xmippdocsym('Iter_1_current_angles.doc','sym_doc.doc','c2','upright'); %without a phys copy of the particle 
%  
%  tom_xmippdocsym('Iter_25_current_angles.doc','sym_doc.doc','c2','normal','/test/target_fold','lns');
%  %with a duplication of the particle using a symbolic link
%
%REFERENCES
%
%SEE ALSO
%   tom_av2_em_classify3d
%
%   created by FB 08/09/09
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

if (nargin<4)
    flag='normal';
end;

if (nargin<5)
    sym_parts_path='';
end;

if (nargin<6)
    cp_flag='lns';
end;

if (nargin<7)
    htl_filename='';
end;

if (nargin<8)
    dir='high2low';
end;

if (nargin<9)
    sym_parts_path_htl='';
end;


if (strcmp(cp_flag,'lns'))
    dpl_pref='ln -s ';
else
    dpl_pref='cp ';
end;


fprintf('%s ', 'reading current ang doc file...');
st=tom_xmippdocread(input_doc);
fprintf('%s \n', ['...done! ' ]); 


s_num=str2double(sym(2));

imcre=360./s_num;



num_of_entries=size(st,1);

if (strcmp(flag,'upright'))
    
    disp('transforming org doc with 180 90 270 zyz conv!!');
    for ii=1:num_of_entries
        new_ang=tom_sum_rotation_zyz([180 90 270;st(ii).rot st(ii).tilt st(ii).psi]);
        st(ii).rot=new_ang(1);
        st(ii).tilt=new_ang(2);
        st(ii).psi=new_ang(3);
    end;
end;    

if (isempty(htl_filename)==0)
    htl=importdata(htl_filename);
    unix(['cp ' htl_filename ' ' htl_filename '.ext']);
    fp_ht=fopen([htl_filename '.ext'],'a');
    if (strcmp(dir,'high2low'))
        htl_p_idx=htl.data(:,1);
    else
        htl_p_idx=htl.data(:,2);
    end;
end;


max_ref=max([st(:).ref]);
tmp_ang_cl=zeros(max_ref,2);
for i=1:length(st)
    tmp_ang_cl(st(i).ref,:)=[st(i).rot st(i).tilt];
end;

ref_count=max_ref;

htl_error_count=0;

st_new=st;
tic;
for i=1:s_num-1
    tmp=st;
    for ii=1:num_of_entries
        %tmp(ii).rot=tmp(ii).rot+imcre;
        tmp(ii).rot=tmp(ii).rot+180;
        
        if (tmp(ii).rot>180)
            tmp(ii).rot=tmp(ii).rot-360;
        end;
        new_ang_tmp=[tmp(ii).rot tmp(ii).tilt];
        
        idx=find(sum(tmp_ang_cl==repmat(new_ang_tmp,length(tmp_ang_cl),1),2)==2);
        
        
        if (isempty(idx))
            ref_count=ref_count+1;
            tmp(ii).ref=ref_count;
            tmp_ang_cl(tmp(ii).ref,:)=[tmp(ii).rot tmp(ii).tilt];
        else
            tmp(ii).ref=idx(1);
        end;
       
        
        if (isempty(sym_parts_path)==0)
            folder_num=floor((ii-1)./10000)+1;
            if (mod(ii,10000) == 1 || ii==1 )
                warning off; mkdir([sym_parts_path '/fold_'  num2str(folder_num)]); warning on;
            end;
            [aa bb cc]=fileparts(tmp(ii).name);
            
            new_path=[sym_parts_path '/fold_'  num2str(folder_num) '/sym'  bb(1:max(strfind(bb,'_')))  num2str(tmp(ii).part_idx+10000000) '.spi'];
            str=[dpl_pref st(ii).name ' ' new_path ];
            unix(str);
            tmp(ii).name=new_path;
        end;
        
        if (isempty(htl_filename)==0)
            try
                folder_num=floor((ii-1)./10000)+1;
                if (mod(ii,10000) == 1 || ii==1 )
                    warning off; mkdir([sym_parts_path_htl '/fold_'  num2str(folder_num)]); warning on;
                end;
                idx=find([st(ii).part_idx==htl_p_idx]);
                new_name1=[sym_parts_path_htl '/fold_'  num2str(folder_num) '/syma_' num2str(htl.data(idx,1)+1000000) '.spi'];
                new_name2=[sym_parts_path_htl '/fold_'  num2str(folder_num) '/symb_' num2str(htl.data(idx,2)+1000000) '.spi'];
                org_name1=htl.textdata{idx,1};
                org_name2=htl.textdata{idx,2};
                
                str1=['cp  ' org_name1 ' ' new_name1 ];
                unix(str1);
                str2=['cp  ' org_name2 ' ' new_name2 ];
               % unix(str2);
                fprintf(fp_ht,'%s %s %s %s\n',new_name1,new_name2, num2str(htl.data(idx,1)+1000000), num2str(htl.data(idx,2)+1000000));
         catch ME
                htl_error_count=htl_error_count+1;
                disp(['htl error in ' num2str(st(ii).part_idx) ]);
                disp([ME.message]);
            end;
        end;
        
        
        if (mod(ii,2000)==0)
            disp(num2str(ii));
            toc;
            tic;
        end;
        
    end;
    st_new=cat(1,st_new,tmp);
end;


if (isempty(htl_filename)==0)
    fclose(fp_ht);
    disp(['htl Errors: ' num2str(htl_error_count)]);
end;

st_new(1).part_idx_unique=0;
fprintf('%s ', 'writing sym doc file...');
tom_xmippdocwrite(output_doc,st_new);
fprintf('%s \n', ['...done! ' ]); 






% [a b c]=fileparts(tmp(1).name);
% tmp_pos=strfind(a,'/');
% fold_name=a(tmp(length(tmp_pos)-3):tmp(length(tmp_pos)-2));







































