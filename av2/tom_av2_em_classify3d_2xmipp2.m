function tom_av2_em_classify3d_2xmipp2(fpart_st,outputfold,class_nr,f_output_sel,f_output_doc,find_what,replace_with)
%TOM_AV2_EM_CLASSIFY3D_2XMIPP2 creates xmipp doc and sel
%   
%
%  tom_av2_em_classify3d_2xmipp2(cl_struct,class_nr,f_output_sel,f_output_doc,find_what,replace_with)
%
%  TOM_AV2_EM_CLASSIFY3D_2XMIPP2 creates xmipp doc and sel (converter
%  function for tom_av2_em_classify3d) for angular refinement with xmipp
%  
%
%PARAMETERS
%
%  INPUT
%   fpart_st                      particle structure (normally runXX/data/part_st.mat)
%   outputfold                    (opt) default ./  
%   class_nr                      (opt) class number (default 'all')
%   out_sel                       (opt) output base for sel (default 'parts_')                     
%   out_doc                       (opt) output base for doc (default 'parts_')
%   find_what                     (opt) string 2 find 
%   replace_with                  (opt) string 2 be replaced                   
%  OUTPUT
%
%EXAMPLE
%  
%  %builds sel and doc for all classes and dumps it 2 xmipp_run1
%  tom_av2_em_classify3d_2xmipp('part_st.mat','xmipp_run1');
%
%  %builds sel and doc for class nr 1,3 and dumps it under xmipp/cl1.selxmipp/cl1.doc  
%   and xmipp/cl1.sel xmipp/cl1.doc 
%  tom_av2_em_classify3d_2xmipp('part_st.mat','xmipp',[1 3],'cl1sel','cl1doc');
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


if (nargin < 2)
    outputfold='./';
end;

if (nargin < 3)
    class_nr='all';
end;

if (nargin < 4 )
    f_output_sel='parts_';
end;

if (nargin < 5 )
    f_output_doc='parts_';
end;

if (nargin<6)
    find_what='';
end;

if (nargin<7)
    replace_with='';
end;

if (strcmp(outputfold,'./')==0)
    warning off;
    mkdir(outputfold);
    warning on;
end;

%load master
part_st_master=load(fpart_st);
part_st_master=part_st_master.st;

xmipp_doc=part_st_master.doc_name;

if (strcmp(class_nr,'all'))
    class_nr=1:part_st_master.num_of_classes;
end;


disp(['Extracting class nr: ' num2str(class_nr) ' from ' num2str(part_st_master.num_of_classes) ' classes']);


for i=1:length(class_nr)
    all_path_out_doc{i}=[outputfold '/' f_output_doc num2str(class_nr(i)) '.doc'];
    all_path_out_sel{i}=[outputfold '/' f_output_sel num2str(class_nr(i)) '.sel'];
end;

fprintf('%s',['reading ' xmipp_doc]);
doc_all=tom_xmippdocread(xmipp_doc);
fprintf('%s \n',' ...done!');


disp(['Found ' num2str(length(doc_all)) ' particles in '  num2str(length(part_st_master.split_part_mat)) ' chunk(s)']);
disp(' ');

all_doc='';
for ii=1:length(class_nr)
     ch_count=1;
     disp(['Processing Class nr ' num2str(class_nr(ii))]); 
    for i_chunk=1:length(part_st_master.split_part_mat)
        part_st=load(part_st_master.split_part_mat{i_chunk});
        part_st=part_st.st;
        disp(['  Processing Chunk nr ' num2str(i_chunk) ' (iter: ' num2str(size(part_st.class,1)) ')' ]);
        class_count=0;
        doc=doc_all(ch_count:(ch_count+length(part_st.ref_nr)-1));
        new_doc=doc;
        classes=part_st.class(size(part_st.class,1),:);
        for i=1:size(part_st.class,2)
            if (isempty(find(class_nr(ii)==classes(i),1))==0)
                if isempty(f_output_doc)==0
                    class_count=class_count+1;
                    new_doc(class_count)=doc(i);
                end;
            end;
        end;
        all_doc{i_chunk,ii}=new_doc(1:class_count);
        ch_count=ch_count+length(part_st.ref_nr);
    end;
end;
clear('new_doc');

disp('Merge Chunks');
for i=1:length(class_nr)
    doc_tmp='';
    for i_chunk=1:length(part_st_master.split_part_mat)
        doc_tmp=cat(1,doc_tmp,all_doc{i_chunk,i});
    end;
    doc_tmp(1).header=doc_all(1).header;
    doc_tmp(1).part_idx_unique=1;
    disp(['  writing ' all_path_out_doc{i} ' (' num2str(length(doc_tmp)) ' parts)' ]);
    tom_xmippdocwrite(all_path_out_doc{i},doc_tmp);
    disp(['  writing ' all_path_out_sel{i}  ' (' num2str(length(doc_tmp)) ' parts)']);
    fid=fopen(all_path_out_sel{i},'wt');
    for ii=1:length(doc_tmp)
        fprintf(fid,'%s 1\n',doc_tmp(ii).name);
    end;
    fclose(fid);
end;
disp('done!');


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
