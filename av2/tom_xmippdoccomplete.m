function comp_flag=tom_xmippdoccomplete(doc_name,new_doc_name,sym_parts_path,flag)
%  tom_xmippdoccomplete completes a doc file from tom_xmippdocsym
%   for xmipp refinement
%  
%    tom_xmippdocsym(input_doc,output_doc,sym)
%  
%    TOM_XMIPPDOCSYM is needed due to the fact that xmipp cannot handle the same particle with different angles 
%                    
%                    
%  PARAMETERS
%  
%    INPUT
%     doc_name               *.doc filename 
%     new_doc_name          name for output doc name 
%     sym_parts_path        path for the sym particels
%     flag                  (lns) or cp
%    
%  
%    OUTPUT
%     comp_flag              flag for incomlete doc 
%                            0=complete 1=incomplete 
%  EXAMPLE
%       
%    tom_xmippdoccomplete('model_1.doc','test_compl.doc','/fs/scratch/fbeck/all_pombe/try_rpn13/test_refine/out_test/model_1/sym_parts/');
%  
%  NOTE
%  
%      
%  
%  REFERENCES
%  
%  SEE ALSO
%     tom_av2_em_classify3d
%  
%     created by FB 08/09/09
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

if (nargin<4)
    flag='ln-s';
end;

comp_flag=0;

%read doc file

disp(['reading ' doc_name]);
doc=tom_xmippdocread(doc_name);
disp('done!');

tmp_nums=[doc(:).part_idx];

[tmp_sort idx_sort]=sort(tmp_nums);
idx=find(diff(tmp_sort)==0);

%allocate memory
new_doc=doc;

%remove identical names
disp(['comleting doc file ' num2str(length(idx)) ' of ' num2str(length(doc)) ' needed' ]);

tic;
for i=1:length(idx)
    old_name=doc(idx_sort(idx(i)+1)).name;
    [a b c]=fileparts(old_name);
    [tok rest]=strtok(b,'_');
    fold_num=floor((i-1)./40000)+1;
    new_name=[sym_parts_path '/fold_' num2str(fold_num) '/' tok 'Sym' rest c];
    
    if (mod(i,40000) == 1 || i==1 || exist([sym_parts_path '/fold_' num2str(fold_num)],'dir')==0 )
        warning off; mkdir([sym_parts_path '/fold_' num2str(fold_num)]); warning on;
    end;
    
    call=[old_name ' ' new_name];
    if (strcmp(flag,'ln-s'))
        call=['ln -s ' old_name ' ' new_name];
    else
        call=['cp ' old_name ' ' new_name];
    end;
    
    [error status]=unix(call);
    if (error==1 && isempty(strfind(status,': File exists'))==1 )
        disp(['error processing: ' doc(idx_sort(idx(i)+1)).name]);
        error(status);
    end;
    new_doc(idx_sort(idx(i)+1)).name=new_name;
    
    if (mod(i,1000)==0)
        toc;
        disp([num2str(i) ' particles processed!']);
        tic;
    end;
    
end;

if (length(idx)>0)
    comp_flag=1;
end;

disp(['writing ' new_doc_name]);
tom_xmippdocwrite(new_doc_name,new_doc);
disp('done');



disp(' ');










