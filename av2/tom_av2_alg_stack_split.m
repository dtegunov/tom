function tom_av2_alg_stack_split(in_stack_alg,in_struct_alg,out_stack_alg,out_struct_alg,out_filt_doc,class_vect,itr)
% TOM_AV2_ALG_STACK_SPLIT splits a aligned and classified stack according 2 class vect
% 
% 
% tom_av2_alg_stack_split(in_stack_alg,in_struct_alg,out_stack_alg,out_struct_alg,class_vect)
%  
% 
% PARAMETERS
% 
% INPUT
% in_stack_alg         input stack name 
% in_struct_alg        input struct (from tom_av2_em_classify3d) name
% out_stack_alg        output stack name
% out_struct_alg       output struct name
% out_filt_doc         output doc name 
% class_vect           vector containing the selected classes
% itr                  (last iter) iteration 2 use                                                 
% 
% 
% OUTPUT
% 
% EXAMPLE
% 
% tom_av2_alg_stack_split('st_cutf2.em','rclLid5/part_st.mat','st_cl1_2.em','st_cl1_2.mat','doc_cl1_2.doc',[1 2]);
% 
% 
% REFERENCES
% 
% SEE ALSO
% tom_av2_em_classify3d,tom_av2_xmipp_align_stack,tom_av2_em_classify3d_2xmipp
% 
% created by FB 08/09/09
% 
% Nickell et al., 'TOM software toolbox: acquisition and analysis for electron tomography',
% Journal of Structural Biology, 149 (2005), 227-234.
% 
% Copyright (c) 2004-2007
% TOM toolbox for Electron Tomography
% Max-Planck-Institute of Biochemistry
% Dept. Molecular Structural Biology
% 82152 Martinsried, Germany
% http://www.biochem.mpg.de/tom


max_chunk_size=250000;

%read struct
load(in_struct_alg);

if (nargin < 7)
    itr=size(part_st.class,1);
end;

disp(['itr: ' num2str(itr) ' is used']);


%read stack
if (isfield(part_st,'split_part_stack'))
    disp('Reading sub stacks: ');
    part_stack=[];
    for i=1:length(part_st.split_part_stack)
       disp(['   ' part_st.split_part_stack{i}]);
       tmp=tom_emreadc(part_st.split_part_stack{i});
       tmp=single(tmp.Value);
       part_stack=cat(3,part_stack,tmp);
       clear('tmp');
    end;
else
    part_stack=tom_emreadc(in_stack_alg);
    part_stack=single(part_stack.Value);
end;


%generate used part idx
idx=find(ismember(part_st.class(itr,:),class_vect));

disp(['class ' num2str(class_vect) ' : ' num2str(length(idx)) ' of ' num2str(size(part_stack,3)) ]);


%filter stack
part_stack=single(part_stack(:,:,idx));
if (size(part_stack,3) > max_chunk_size )
    disp(['lenght stack: ' num2str(size(part_stack,3)) ' > ' num2str(max_chunk_size) ' ==> splitting stack'  ])
    [f_base f_name f_ext]=fileparts(out_stack_alg);
    packages=tom_calc_packages(ceil(size(part_stack,3)./max_chunk_size),size(part_stack,3));
    for ip=1:size(packages,1)
        if (isempty(f_base))
            f_base=pwd();
        end;
        tmp_name=[f_base '/' f_name '_' num2str(ip) f_ext];
        split_part_stack{ip}=tmp_name;
        disp(['   stack alg ==> ' num2str(packages(ip,1)) ' to ' num2str(packages(ip,2)) ' ==> ' tmp_name] );
        tom_emwritec(tmp_name,part_stack(:,:,packages(ip,1): packages(ip,2)),'standard','single');
    end;
else
    tom_emwritec(out_stack_alg,part_stack,'standard','single');
end;

clear('part_stack');


%filter doc
doc=tom_xmippdocread(part_st.doc_name);
tmp_head=doc(1).header;
tmp_unique=doc(1).part_idx_unique;
doc=doc(idx);
doc(1).header=tmp_head;
doc(1).part_idx_unique=tmp_unique;
tom_xmippdocwrite(out_filt_doc,doc);


%filter struct
new.ref_nr=part_st.ref_nr(idx);
new.euler_ang_zxz=part_st.euler_ang_zxz(idx);
new.euler_ang_zyz=part_st.euler_ang_zyz(idx);
new.error_idx=part_st.error_idx(idx);
new.part_names=part_st.part_names(idx);
new.doc_name=out_filt_doc;
new.error_idx=part_st.error_idx(idx);
new.euler_ang_zxz_proj=part_st.euler_ang_zxz_proj;
new.euler_ang_zyz_proj=part_st.euler_ang_zyz_proj;
new.pre_filter_k=part_st.pre_filter_k;
new.split_part_stack=split_part_stack;
part_st=new;
for i=1:length(part_st.part_names)
    part_st.part_names{i}=strrep(part_st.part_names{i},'.//fs/pool/','/fs/pool/');
end;
part_st.org_part_st=in_struct_alg;

save(out_struct_alg,'part_st');
clear('part_st');












