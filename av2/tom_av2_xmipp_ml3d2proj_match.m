function tom_av2_xmipp_ml3d2proj_match(basep_iter_num,target_folder,htl_filename,all_flag,find_what,replace_with,proj_match_scrt,f_ctfdat,f_defocus_doc,class_nr)
%TOM_AV2_XMIPP_ML3D2PROJ_MATCH prepares files for xmipp-projmatch
%
%   tom_av2_xmipp_ml3d2proj_match(basep_iter_num,target_folder,htl_filename,all_flag,find_what,replace_with,proj_match_scrt)
%
%  TOM_AV2_XMIPP_ML3D2PROJ_MATCH eprepares files for xmipp-projmatch
%  after ml3d classification
%  
%  
%
%PARAMETERS
%
%  INPUT
%   basep_iter_num     basep with iteration number  (model_it000002)
%   target_folder      folder files should be transferd use abs path 
%   htl_filename       (opt) high2low filename
%   all_flag           (opt) get high and low (default 0)
%   find_what          (opt) string 2 find 
%   replace_with       (opt) string 2 be replaced                     
%   proj_match_scrt    (opt) xmipp proj match script 
%   f_ctfdat           (opt) filename of xmipp .ctfdat file
%   f_defocus_doc      (opt) filename of xmipp doc-file 2 split in defocus groups
%   class_nr           (opt) vectro with class numbers to process (default all)
%
%EXAMPLE
%    
% %%% WORKING WITH TOM_AV2_EM_CLASSIFY3D %%%%%%%%
%  
%
% tom_av2_xmipp_ml3d2proj_match('part_st.mat','/fs/pool/pool-nickell2/scratch/te_d/out','',0,'','','xmipp_protocol_projmatch.py');
% 
%
%
%
% %%% WORKING WITH XMIPP ML3D %%%%%%%%
%
% tom_av2_xmipp_ml3d2proj_match('out_ml3d/model_it000011','/fs/pool/pool-nickell2/scratch/classify_from_blue_gene/ml3d_split_two_5/rCtfg','11__20_high2low.htl',0,'/u/fbeck/BlueGene/data_bohn/','/fs/scratch/bohn/','xmipp_protocol_projmatch.py','11__20_corrf_1st_128_ctfDat.txt','ctf_groups_split.doc');
% %all you can eat example which performs a particleExtraction, path_addaption, high2low, ctf-groups and a polymorph script-generation 
%
%
% tom_av2_xmipp_ml3d2proj_match('out_ml3d/model_it000011','/fs/pool/pn/refine','11__20_high2low.htl',0,'/u/fbeck/BlueGene/data_bohn/','/fs/scratch/bohn/','xmipp_protocol_projmatch.py');
% % example which performs a particleExtraction, path_addaption, high2low addaption and a polymorph script-generation 
%
%
% tom_av2_xmipp_ml3d2proj_match('out_ml3d/model_it000011','/fs/pool/pn/refine','11__20_high2low.htl',1,'/u/fbeck/BlueGene/data_bohn/','/fs/scratch/bohn/');  
% %example without script generation !
%
%  %get data from BlueGene!
%
%  %from unix shell sync with:
%
%  %log in  ...to open port 6001
%  ssh -L localhost:6001:genius1:22 fbeck@gate.rzg.mpg.de 
%  ssh genius
%  %sync it  (...open new shell)
%  rsync  -avz -e 'ssh -p '"6001"
%  "fbeck@localhost:./data_bohn/log/ml3d_all_high_64_sort_rubbish/out" "./" --exclude "*.xmp" --exclude "*.proj"
%   
%  %volums only
%  rsync  -avz -e 'ssh -p '"6001" "fbeck@localhost:./data_bohn/log_270k/ml3d_all_high_64_sort_rubbish/out/model_it000010_vol0000??.vol" "."
%
%
%NOTE
% 
% !!! path for target_folder shoul be absolute !!!!
% !!!abs path for target_folder shoul be loner than 70 characters (>70 xmipp crashes with seg. fault in Reconstruction Module) !!!! 
% % 
% % proj_match_scrt rules (start volume,doc file name ,ctf_dat and ctf-group files have to have the following names !)
% 
%       *************************************************  
%       *  start volume    ==>   start.spi              *
%       *  doc file name   ==>   my_doc.doc (expert o.) *      
%       *  ctf_dat         ==>   all_images.ctfdat      *
%       *  defoucus split  ==>   split_ctf.doc          *
%       *************************************************
%
%REFERENCES
%
%SEE ALSO
%   tom_av2_xmipp_ml3d_part_per_class,tom_av2_em_classify3d_2xmipp,tom_av2_xmipp_ml3d2proj_match_scrt_adapt,tom_xmippdocread,
%
%   created by fb 
%
%   Nickell et al., 'TOM software toolbox: acquisition and analysis for electron ponography',
%   Journal of Structural Biology, 149 (2005), 227-234.
%
%   Copyright (c) 2004-2007
%   TOM toolbox for Electron Tomography
%   Max-Planck-Institute of Biochemistry
%   Dept. Molecular Structural Biology
%   82152 Martinsried, Germany
%   http://www.biochem.mpg.de/tom


if (nargin<3)
    htl_filename='';
end;

if (nargin<5)
    all_flag=0;
end;

if (nargin<5 || nargin<6)
    find_what='';
end;

if (nargin<7)
    proj_match_scrt='';
end;

if (nargin<8)
    f_ctfdat='';
end;

if (nargin<9)
    f_defocus_doc='';
end;

if (nargin<10)
    class_nr='';
end



beck_flag=1;
try
    load(basep_iter_num);
catch ME
    beck_flag=0;
end;

if (beck_flag==0)
    %retrive information from ml3d output folder
    [base_path base_name]=fileparts(basep_iter_num);
    all_volumes=dir([basep_iter_num '*.vol']);
    doc=dir([basep_iter_num '*.doc']);
    all_sels=dir([basep_iter_num '_*.sel']);
    num_of_classes=length(all_volumes);
    [rr ar]=strtok(base_name,'it');
    ar=strrep(ar,'it','');
    max_iter=str2num(ar);
else
    %retrive information from ml3d output mat file
    max_iter=size(part_st.class,1);
    %[base_path base_name]=fileparts(part_st.outputdir);
    base_path=part_st.outputdir;
    base_path=[base_path '/models/'];
    all_volumes=dir([part_st.outputdir '/models/model_iter_' num2str(max_iter)  '_model*.em']);
    if (isempty(all_volumes))
        max_iter=max_iter-1;
        all_volumes=dir([part_st.outputdir '/models/model_iter_' num2str(max_iter)  '_model*.em']);
    end;
    num_of_classes=length(all_volumes);
end;


if (isempty(htl_filename)==0 )
    disp(['reading ' htl_filename]);
    try
        in_htl=importdata(htl_filename);
    catch ME
        disp(['Cannot open ' htl_filename]);
        error(ME.message);
    end;
    im_tmp=tom_spiderread(in_htl.textdata{1});
    vol_size=[size(im_tmp.Value) size(im_tmp.Value,1)];
    clear('low_path'); clear('im_tmp');
else
    if (beck_flag==0)
    else
        doc=tom_xmippdocread(part_st.doc_name);
        im_tmp=tom_spiderread(doc(1).name);
        vol_size=[size(im_tmp.Value) size(im_tmp.Value,1)];
        clear('im_tmp');
    end;
end;


if (isempty(class_nr))
    class_nr=1:num_of_classes;
end;

if (beck_flag==1)
    disp(['3D EM (TOM) Cl with ' num2str(num_of_classes) ' classes and ' num2str(max_iter) ' iterations found' ]);
else
    disp(['ML3d Cl (Xmipp) with ' num2str(num_of_classes) ' classes and ' num2str(max_iter) ' iterations found' ]);
end;
disp(' ');

if (exist(target_folder,'dir')==0)
    mkdir(target_folder);
end;

fid=fopen([target_folder '/st.sh'],'w');

if (num_of_classes==0)
    error(['no volumes found check ' basep_iter_num]);
end;

%main loop over all classes
for i=class_nr
    
    disp(' '); disp(['Processing Class nr: ' num2str(i)]); disp(' ');
    out_basepath=[target_folder '/model_' num2str(i) ];
    warning off;
    mkdir(out_basepath);
    mkdir([out_basepath '/org']);
    warning on;
    
    in_start_volume=[base_path '/' all_volumes(i).name];
    out_start_volume=[out_basepath '/start.spi'];
    
    out_sel=[out_basepath '/org/model_' num2str(i) '.sel'];
    out_doc=[out_basepath '/org/model_' num2str(i) '.doc'];
    
    disp(' '); disp(['Starting Particle extraction... ' ]); disp(' ');
    
    if (beck_flag==0)
        % extract from xmipp ml3d folder
        in_sel=[base_path '/' all_sels(i).name];
        in_doc=[base_path '/' doc(1).name];
        if (isempty(find_what))
            tom_av2_xmipp_ml3d_part_per_class(in_doc,in_sel,out_sel,out_doc);
        else
            tom_av2_xmipp_ml3d_part_per_class(in_doc,in_sel,out_sel,out_doc,find_what,replace_with);
        end;
    else
        % extract from tom 3d em output folder
        if (isempty(find_what))
            tom_av2_em_classify3d_2xmipp(basep_iter_num,i,out_sel,out_doc);
        else
            tom_av2_em_classify3d_2xmipp(basep_iter_num,i,out_sel,out_doc,find_what,replace_with);
        end;
    end;
    
    
    if (isempty(find_what)==0)
        [base_s name_s ext_s]=fileparts(out_sel);
        [base_d name_d ext_d]=fileparts(out_doc);
        name_sel=[base_s '/' name_s '_path' ext_s];
        name_doc=[base_d '/' name_d '_path' ext_d];
    end;
    
    if (isempty(htl_filename) )
        % no htl file so just copy corresponding files from org !!
        if (isempty(find_what))
            unix(['cp ' out_sel ' ' out_basepath '/model_' num2str(i) '.sel']);
            unix(['cp ' out_doc ' ' out_basepath '/model_' num2str(i) '.doc']);
        else
            unix(['cp ' name_sel ' ' out_basepath '/model_' num2str(i) '.sel']);
            unix(['cp ' name_doc ' ' out_basepath '/model_' num2str(i) '.doc']);
        end;
        if (tom_isemfile(in_start_volume)==0)
            vol_tmp=tom_spiderread(in_start_volume);
        else
            vol_tmp=tom_emread(in_start_volume);
        end;
        vol_tmp=tom_rescale3d(vol_tmp.Value,vol_size);
        tom_spiderwrite(out_start_volume,vol_tmp);
        
        
    else
        % feed tom_av2_high2low_transform2 ... to get corresponding low files !
        disp(' '); disp('Starting high2low transform... ' ); disp(' ');
        if (isempty(find_what))
            tom_av2_high2low_transform2(out_sel,out_doc,htl_filename,[out_basepath '/model_' num2str(i) '.sel'],[out_basepath '/model_' num2str(i) '.doc'],'high2low',all_flag);
        else
            tom_av2_high2low_transform2(name_sel,name_doc,htl_filename,[out_basepath '/model_' num2str(i) '.sel'],[out_basepath '/model_' num2str(i) '.doc'],'high2low',all_flag);
        end;
        if (tom_isemfile(in_start_volume)==0)
            vol_tmp=tom_spiderread(in_start_volume);
        else
            vol_tmp=tom_emread(in_start_volume);
        end;
        vol_tmp=tom_rescale3d(vol_tmp.Value,vol_size);
        tom_spiderwrite(out_start_volume,vol_tmp);
    end;
    
    %complete doc file
    comp_flag=tom_xmippdoccomplete([out_basepath '/model_' num2str(i) '.doc'],[out_basepath '/model_' num2str(i) '.doc'],[out_basepath '/tmp_parts/']);
    if (comp_flag==1)
        tom_xmippdoc2sel([out_basepath '/model_' num2str(i) '.doc'],[out_basepath '/model_' num2str(i) '.sel']);
        disp('Rotateing doc file 0 90 90 in zxz');
        tom_xmippdocrotate([out_basepath '/model_' num2str(i) '.doc'],[out_basepath '/model_' num2str(i) '.doc'],[0 90 90]);
        tom_spiderwrite(out_start_volume,tom_rotate(vol_tmp,[-90 0 -90]));
        disp('done');
    end;
    
    
    if (isempty(proj_match_scrt)==0)
        try
            tom_av2_xmipp_ml3d2proj_match_scrt_adapt(proj_match_scrt,target_folder,i);
            [ps_p ps_n ps_e]=fileparts(proj_match_scrt);
            new_proj_path=out_basepath;
            fprintf(fid,['cd ' new_proj_path '\n']);
            fprintf(fid,[new_proj_path '/' ps_n ps_e '\n']);
            fprintf(fid,['chmod -R ugo+rwx' ' ' new_proj_path '&\n']);
        catch ME
            disp(ME.message);
            disp('Error adapting xmipp script ...skipping (try with tom_av2_xmipp_ml3d2proj_match_scrt_adapt again!!)');
        end;
    end;
    
    %copy ctfdat and defocus.doc
    if (isempty(f_ctfdat)==0)
        disp(['copy ' f_ctfdat ' ==> ' out_basepath '/all_images.ctfdat']);
        unix(['cp ' f_ctfdat ' ' out_basepath '/all_images.ctfdat']);
    end;
    
    if (isempty(f_defocus_doc)==0)
        disp(['copy ' f_defocus_doc ' ==> ' out_basepath '/split_ctf.doc']);
        unix(['cp ' f_defocus_doc ' ' out_basepath '/split_ctf.doc']);
    end;
    
end;


fclose(fid);
unix(['chmod ugo+rwx ' target_folder '/st.sh']);
copyfile([target_folder '/st.sh'],[target_folder '/start.sh'])
unix(['rm ' target_folder '/st.sh']);
















