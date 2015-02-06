function tom_av2_xmipp_ml3d_part_per_class(doc_file,class_sel,f_output_sel,f_output_doc,find_what,replace_with)
%TOM_AV2_XMIPP_ML3D_PART_PER_CLASS extracts particles according to a class sel
%
%   tom_av2_xmipp_ml3d_part_per_class(doc_file,class_sel,f_output_sel,f_output_doc,find_what,replace_with)
%
%  tom_av2_xmipp_ml3d_part_per_class extracts particles according to a class sel
%  and creates a sel and a doc file ...the root path can also be changed 
%  
%  
%
%PARAMETERS
%
%  INPUT
%   doc_file         filename of the input doc
%   class_sel        filename of the class sel
%   f_output_sel     output sel filename
%   f_output_doc     output doc-file name
%   find_what        (opt) string 2 find 
%   replace_with     (opt) string 2 be replaced                     
%
%EXAMPLE
%     
%    
% tom_av2_xmipp_ml3d_part_per_class('model_it000016_small.doc','model_it000016_vol000001.sel','out_cl1_fb.sel','out_cl1_fb.doc','/u/fbeck/BlueGene/data_bohn/','/fs/scratch/bohn/');
% %all you can eat example creates a sel and doc and addapts the path 
%
% tom_av2_xmipp_ml3d_part_per_class('model_it000016_small.doc','model_it000016_vol000001.sel','out_cl1_fb.sel','out_cl1_fb.doc');
% %minimal example
%
%
%REFERENCES
%
%SEE ALSO
%  tom_xmippdocread
%
%   created by fb ...ole !!
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

if (nargin < 5)
    find_what='';
end;

if (nargin < 4)
    f_output_doc='';
end;


if (strcmp(class_sel,'all')==0)
    
    if (isempty(f_output_doc)==0)
        unix(['head -1 ' doc_file ' > ' f_output_doc]);
    end;
    
    disp('reading cell...');
    sel=importdata(class_sel);
    
    tmp_st=sel.textdata{1};
    num=findstr(sel.textdata{1},'ref');
    tmp_st=tmp_st(num+3:end);
    tmp_st=strrep(tmp_st,'.xmp','');
    start=str2num(tmp_st);
    
    tmp_st=sel.textdata{end};
    num=findstr(sel.textdata{1},'ref');
    tmp_st=tmp_st(num+3:end);
    tmp_st=strrep(tmp_st,'.xmp','');
    stop=str2num(tmp_st);
    
    disp('reading doc...');
    doc=tom_xmippdocread(doc_file);
    
    disp(' ');
    
    disp(['extracting from ' doc_file  ' ==> ' f_output_sel]);
    fp=fopen(f_output_sel,'wt');
    zz=1;
    tic;
    used_idx=zeros(length(doc),1);
    zz_u_idx=0;
    for i=1:length(doc)
        if (doc(i).ref>=start && doc(i).ref<=stop)
            fprintf(fp,[doc(i).name ' 1 \n']);
            zz=zz+1;
            if (isempty(f_output_doc)==0)
                %call=['grep -A 1 "' doc(i).name '" '  doc_file ' >> ' f_output_doc];
                %unix(call);
                zz_u_idx=zz_u_idx+1;
                used_idx(zz_u_idx)=i;
            end;
        end;
        if (mod(i,10000)==0)
            disp([num2str(i) ' of ' num2str(length(doc)) ' done!' ]);
            toc;
            tic;
        end;
    end;
    fclose(fp);
    used_idx=used_idx(1:zz_u_idx);
    
    new_doc=doc(used_idx);
    new_doc(1).part_idx_unique=doc.part_idx_unique;
    new_doc(1).header=doc(1).header;
    tom_xmippdocwrite(f_output_doc,new_doc);
    
    
    disp(' ');
    disp([num2str(zz-1) ' of ' num2str(length(doc)) ' in req class']);
    
    disp(' ');
    
    
    if (isempty(find_what)==0)
        
        [a b c]=fileparts(f_output_sel);
        if (isempty(a))
            a='.';
        end;
        disp(['adapting path: ' f_output_sel ' > ' find_what ' ==> ' replace_with ' > ' a '/' b '_path' c]);
        call=['awk ''{gsub("'  find_what '","' replace_with '"); print }'' ' f_output_sel ' > ' a '/' b '_path' c];
        unix(call);
        if (isempty(f_output_doc)==0)
            [a b c]=fileparts(f_output_doc);
            if (isempty(a))
                a='.';
            end;
            disp(['adapting path: ' f_output_doc ' > ' find_what ' ==> ' replace_with ' > ' a '/'  b '_path' c]);
            call=['awk ''{gsub("'  find_what '","' replace_with '"); print }'' ' f_output_doc ' > '  a '/' b '_path' c];
            unix(call);
        end;
        
        
    end;
    
    
else
    [a b]=fileparts(doc_file);
    d=dir([a '/' b '_vol*.sel']);
    
    for i=1:length(d)
        all_sels{i}=[a '/' d(i).name];
    end;
    fp=fopen(f_output_sel,'wt');
    disp('reading doc...');
    doc=tom_xmippdocread(doc_file);
    
    for ii=1:length(all_sels)
        class_sel=all_sels{ii};
        
        disp('reading cell...');
        sel=importdata(class_sel);
        
        tmp_st=sel.textdata{1};
        num=findstr(sel.textdata{1},'ref');
        tmp_st=tmp_st(num+3:end);
        tmp_st=strrep(tmp_st,'.xmp','');
        start=str2num(tmp_st);
        
        tmp_st=sel.textdata{end};
        num=findstr(sel.textdata{1},'ref');
        tmp_st=tmp_st(num+3:end);
        tmp_st=strrep(tmp_st,'.xmp','');
        stop=str2num(tmp_st);
        
      
        
        disp(' ');
        
        disp(['extracting from ' doc_file  ' ==> ' f_output_sel]);
       
        zz=1;
        tic;
        used_idx=zeros(length(doc),1);
        zz_u_idx=0;
        for i=1:length(doc)
            if (doc(i).ref>=start &&  doc(i).ref<=stop)
                [a b c]=fileparts(doc(i).name);
                [a d]=strtok(b,'_');
                num=(strrep(d,'_',''));
                
                fprintf(fp,[num ' ' num2str(ii) '\n']);
                zz=zz+1;
                if (isempty(f_output_doc)==0)
                    %call=['grep -A 1 "' doc(i).name '" '  doc_file ' >> ' f_output_doc];
                    %unix(call);
                    zz_u_idx=zz_u_idx+1;
                    used_idx(zz_u_idx)=i;
                end;
            end;
            if (mod(i,10000)==0)
                disp([num2str(i) ' of ' num2str(length(doc)) ' done!' ]);
                toc;
                tic;
            end;
        end;
      
       % used_idx=used_idx(1:zz_u_idx);
        
      
         disp(' ');
        disp([num2str(zz-1) ' of ' num2str(length(doc)) ' in req class']);
        
        disp(' ');
    end;
      fclose(fp);
    
end;




