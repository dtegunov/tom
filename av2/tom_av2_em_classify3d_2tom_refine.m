function tom_av2_em_classify3d_2tom_refine(basep_iter_num,target_folder,htl_filename,all_flag,class_nr)
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
%   class_nr           (opt) class nr for refinement
%
%EXAMPLE
%    
% %%% WORKING WITH TOM_AV2_EM_CLASSIFY3D %%%%%%%%
%  

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

if (nargin<4)
    all_flag=0;
end;

if (nargin<5)
    class_nr='';
end




load(basep_iter_num);




%retrive information from ml3d output mat file
max_iter=size(part_st.class,1);
base_path=part_st.outputdir;
base_path=[base_path '/models/'];
all_volumes=dir([part_st.outputdir '/models/model_iter_' num2str(max_iter)  '_model*.em']);
if (isempty(all_volumes))
    max_iter=max_iter-1;
    all_volumes=dir([part_st.outputdir '/models/model_iter_' num2str(max_iter)  '_model*.em']);
end;
num_of_classes=length(all_volumes);



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
end;
 
if (isempty(class_nr))
    class_nr=1:num_of_classes;
end;

disp(['3D EM (TOM) Cl with ' num2str(num_of_classes) ' classes and ' num2str(max_iter) ' iterations found' ]);

disp(' ');


if (exist(target_folder,'dir')==0)
    mkdir(target_folder);
end;



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
    
    out_sel=[out_basepath '/org/model_' num2str(i) '.sel'];
    out_doc=[out_basepath '/org/model_' num2str(i) '.doc'];
    
    disp(' '); disp(['Starting Particle extraction... ' ]); disp(' ');
    %out1
    tom_av2_em_classify3d_2xmipp(basep_iter_num,i,out_sel,out_doc);
    
    if (isempty(htl_filename) )
        % no htl file so just copy corresponding files from org !!
        unix(['cp ' out_sel ' ' out_basepath '/model_' num2str(i) '.sel']);
        unix(['cp ' out_doc ' ' out_basepath '/model_' num2str(i) '.doc']);
    else
        % feed tom_av2_high2low_transform2 ... to get corresponding low files !
        disp(' '); disp('Starting high2low transform... ' ); disp(' ');
        %out2
           tom_av2_high2low_transform2(out_sel,out_doc,htl_filename,[out_basepath '/model_' num2str(i) '.sel'],[out_basepath '/model_' num2str(i) '.doc'],'high2low',all_flag);
    end;
    
    disp(' '); disp('Starting Stack Alignment... ' ); disp(' ');
    %out3
    tom_av2_xmipp_align_stack([out_basepath '/model_' num2str(i) '.doc'],[out_basepath '/st_' num2str(i) '.em'],[out_basepath '/st_' num2str(i) '.mat'],0);
    disp('done');
    
    disp(' '); disp('Starting Refinement... ' ); disp(' ');
    tom_av2_em_refine3d([out_basepath '/st_' num2str(i) '.em'],[out_basepath '/st_' num2str(i) '.mat'],ones(vol_size),1,3,[out_basepath '/ref'],0,'C1',0,'norm',25,20);
    
    disp('done!');
end;









