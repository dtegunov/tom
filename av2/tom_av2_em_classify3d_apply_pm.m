function tom_av2_em_classify3d_apply_pm(part_st,class,doc_high,doc_low,htl,tmp_data_path,org_part_num,ori_flag,output_dir)
%  tom_xmippdocsym syms a doc file
%     file
%  
%    tom_xmippdocsym(input_doc,output_doc,sym)
%  
%    TOM_XMIPPDOCSYM syms a doc file
%                    
%                    
%  PARAMETERS
%  
%    INPUT
%     part_st               input particle struct 
%     class                 selected class
%     doc_high              doc-file with refined high particles (idx are the same)
%     doc_low               doc-file with refined low particels (idx can differ according 2 htl file)
%     htl_filename          high2low filename (will be extended 2 new particles) 
%     org_part_num          original particle number   
%     ori_flag              orientation falg upright positon ('upright') 
%     output_dir            directory for output
%  
%    OUTPUT
%  
%  EXAMPLE
%       
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


load part_st;

max_iter=size(part_st.class,1);
num_of_cl=size(part_st.class,2);
disp(['tom cl with: ' num2str(num_of_cl) ' particles and ' num2str(max_iter) ' iterations']);

part_idx=find(part_st.class(max_iter,1:146348)==class);

fp=fopen('test.sel','wt');
for i=1:length(part_idx)
    fprintf(fp,'%s 1 \n', part_st.part_names{part_idx(i)});
end;    
fclose(fp);

disp(' ');


function do_match(part_nr,dir_flag)


    
    if (strcmp(dir_flag,'high2low'))
        idx=find(names_h_tmp==in_sel_tmp(i),1);
        if (isempty(idx)==0)
            fprintf(fp,[names_l{idx} ' 1\n']);
            in_count=in_count+1;
        else
            out_count=out_count+1;
        end;
        if (all_flag==1 && isempty(idx)==0 && strcmp(names_h{idx},names_l{idx})==0)
            fprintf(fp,[names_h{idx} ' 1\n']);
        end;
        if isempty(idx)==0
            is_name=in_sel.textdata{i};
            corr_name=names_l{idx};
            corr_name2=names_h{idx}; %attached for both flag
        end;
    end;
    
    if (strcmp(dir_flag,'low2high'))
        idx=find(names_l_tmp==in_sel_tmp(i),1);
        if (isempty(idx)==0)
            fprintf(fp,[names_h{idx} ' 1\n']);
            in_count=in_count+1;
        else
            out_count=out_count+1;
        end;
        if (all_flag==1 && isempty(idx)==0 && strcmp(names_h{idx},names_l{idx})==0)
            fprintf(fp,[names_l{idx} ' 1\n']);
        end;
        if isempty(idx)==0
            is_name=in_sel.textdata{i};
            corr_name=names_h{idx};
            corr_name2=names_l{idx}; %attached for both flag
        end;
    end;
    
    
        
    
   
    



disp(' ');





