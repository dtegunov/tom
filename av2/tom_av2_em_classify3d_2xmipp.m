function  tom_av2_em_classify3d_2xmipp(cl_struct,class_nr,f_output_sel,f_output_doc,find_what,replace_with)
%TOM_AV2_EM_CLASSIFY3D_2XMIPP creates xmipp doc and sel
%   
%
%  tom_av2_em_classify3d_2xmipp(cl_struct,class_nr,f_output_sel,f_output_doc,find_what,replace_with)
%
%  TOM_AV2_EM_CLASSIFY3D_2XMIPP creates xmipp doc and sel (converter
%  function for tom_av2_em_classify3d) for angular refinement with xmipp
%  
%
%PARAMETERS
%
%  INPUT
%   cl_struct                     mat-strcut 
%   class_nr                      (opt) class number (default all classes)
%   out_sel                       (opt) output sel filename (default ./xmipp_doc _cl classnr .sel)                       
%   out_doc                       (opt) output doc-file name (default ./xmipp_doc _cl classnr .sel) 
%   find_what                     (opt) string 2 find 
%   replace_with                  (opt) string 2 be replaced                   
%  OUTPUT
%
%EXAMPLE
%  
%  %builds sel and doc for all classes and dumps it local (./) 
%  tom_av2_em_classify3d_2xmipp('part_st.mat');
%
%  %builds sel and doc for class nr 1 and dumps it under xmipp/cl1.sel xmipp/cl1.doc 
%  tom_av2_em_classify3d_2xmipp('part_st.mat',1,'xmipp/cl1.sel','xmipp/cl1.doc');
%
%
%REFERENCES
%
%SEE ALSO
%   tom_av2_em_classify3d,tom_av2_xmipp_align_stack,tom_av2_xmipp_ml3d2proj_match
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

load(cl_struct);


if (nargin<5)
    find_what='';
end;


if (nargin<6)
    replace_with='';
end;


xmipp_doc=part_st.doc_name;

if (nargin<2)
    f_output_doc='dummy';
    f_output_sel='dummy';
    class_nr=unique(part_st.class(end,:));
end;

[a b c]=fileparts(xmipp_doc);


classes=part_st.class(size(part_st.class,1),:);

if (length(class_nr) > 1)
    for i=1:length(class_nr)
        all_path_out_doc{i}=[ './' b '_cl' num2str(class_nr(i)) '.doc'];
        all_path_out_sel{i}=[ './' b '_cl'  num2str(class_nr(i)) '.sel'];
    end;
else
     for i=1:max(classes)
        all_path_out_doc{i}=f_output_doc;
        all_path_out_sel{i}=f_output_sel;
     end;   
end;

disp(['reading ' xmipp_doc]);
doc=tom_xmippdocread(xmipp_doc);
disp(' done!');

index_names=tom_filenames2index(part_st.part_names);

try
    if (sum([doc(:).part_idx]-index_names')~=0 )
        error('doc or part_st corrupted!');
    end;
catch Me
    disp('No check non unique part idx!');
end;



for ii=1:length(class_nr)
    
    class_count=0;
    fp=fopen(all_path_out_sel{class_nr(ii)},'w');
    new_doc=doc;
    for i=1:size(part_st.class,2)
        if (isempty(find(class_nr(ii)==classes(i)))==0)
            fprintf(fp,[part_st.part_names{i} ' 1\n']);
            if isempty(f_output_doc)==0
                class_count=class_count+1;
                new_doc(class_count)=doc(i);
            end;
        end;
    end;
    fclose(fp);
    new_doc=new_doc(1:class_count);
    new_doc(1).header=doc(1).header;
    new_doc(1).part_idx_unique=1;
    tom_xmippdocwrite(all_path_out_doc{class_nr(ii)},new_doc);
end;



if (isempty(find_what)==0 ) && (length(class_nr) > 1)
    error('path adaption not implemented for more than one class');
    return;
end;


if (isempty(find_what)==0 )
    [a b c]=fileparts(f_output_sel);
     disp(['adapting path: ' f_output_sel ' > ' find_what ' ==> ' replace_with ' > ' a '/' b '_path' c]);
     call=['awk ''{gsub("'  find_what '","' replace_with '"); print }'' ' f_output_sel ' > ' a '/' b '_path' c];
     unix(call);
     if (isempty(f_output_doc)==0)
         [a b c]=fileparts(f_output_doc);
         disp(['adapting path: ' f_output_doc ' > ' find_what ' ==> ' replace_with ' > ' a '/'  b '_path' c]);
         call=['awk ''{gsub("'  find_what '","' replace_with '"); print }'' ' f_output_doc ' > '  a '/' b '_path' c];
         unix(call);   
     end;
end;



%%% VINTAGE CODE %%%%%%%%%%%%%%%%%%%%%%%%%%



% for i=1:length(class_nr)
%    unix(['head -1 ' xmipp_doc ' > '  all_path_out_doc{i}]);
% end;
% 
% 
% for i=1:size(part_st.class,2)
%     if (isempty(find(class_nr==classes(i)))==0)
%         call=['grep -A1 ' part_st.part_names{i} ' ' xmipp_doc ' >> '   all_path_out_doc{classes(i)}];
%         unix(call);
%         call=['grep ' part_st.part_names{i} ' ' xmipp_doc ' >> '  all_path_out_sel{classes(i)}];
%         unix(call);
%         if (mod(i,1000)==0)
%             disp([num2str(i) ' ' ]);
%         end;
%     end;
% end;
% 
% %clean .sel files
% for i=1:length(class_nr)
%     call=['cat '  all_path_out_sel{class_nr(i)}  ' | awk ''{print $2 " 1" }''  > ' fileparts(all_path_out_sel{class_nr(i)}) '/tmp.txt'];
%     unix(call);
%     call=['mv '  fileparts(all_path_out_sel{class_nr(i)}) '/tmp.txt '  all_path_out_sel{class_nr(i)}];
%     unix(call);
%     if (tom_av2_xmipp_check_unique(all_path_out_sel{class_nr(i)})==0)
%         error([base_out num2str(i) '.sel' ' not unique!!']);
%     end;
% end;









