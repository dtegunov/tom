function call_pick()

%% Picking

%execute picker
% -0 picking is skipped use 2 for restart
do_pick=1;                 

%output folder
% -folder for final picklists
outputfolder='/fs/pool/pool-doedel/271206/log/pick/autopick';

%input micrograph path
% -wildcard for input images or sell
img_path='/fs/pool/pool-doedel/271206/high/high_*.em';
%img_path{1}='/fs/pool/pool-doedel/img/high_1.em';


%micrograph pixelsize
%pixelsize of the micrograph in Ang
pixelsize_micrograph=2.7;

%path match lib
% -path to templates 4 matching 
path_match_lib='/fs/pool/pool-bmsan-apps/app_data/tom/autopick/26S/lib_match_150_8.8A';

%pixelsize of match lib
% -pixelsize of the match lib usually encoded in the name 
pixelsize_match_lib=8.8;

%invert image 
% -use 1 for cryo 0 for neg. stain mass is by def white
invert=0;

%max number of particles
% -maximum number of hits per micrograph 
max_num_of_particles=500;

%zero radius
% -expert radius for setting peaks in ccf 2 zero in percent 
% -expert 100 means radius=(template_size./2)
% -expert use values smaller the 100 for non globular objects
zero_rad=100;

%ccc thresh 4 matching
% -threshold for matching use -1 for adjusting
% check tom_av2_filter_align for adjusting
cc_threshold=0.26;

%xcf,soc and psr ratio 
% -expert values for linear combination of ccf functions check: hrabe et al 
xcf=0.85;
soc=0.15;
psr=0;

%do multiref
% -switch for multiref alignment 2 get a better classification 
do_multiref=1;

%path match lib
% -path to multiref struct
mulit_path_lib='/fs/pool/pool-bmsan-apps/app_data/tom/autopick/26S/lib_multiref_g25_b0_8.8A.mat';

% cc-threshold for multiref classification (SELECTIVITY !!)
% threshold 2 discriminate between crap and particle
multi_cc_threshold=0.32; % use ~0.2 to 0.5 for cryo ~0.5 to 0.8 neg stain -1 for adjusting

%use refereces
% -use only_good to make the picker less selective
multi_use_refs='all'; 

%max alignment shift
% -expert maximum allowed shift 4 multiref in percent of (template_size./2)
% -expert 100 means max shift=(template_size./2)
alg_max_sh=8;


%multiref filter type
% -expert type of mulitref filter kernel or none 
multi_filter_type='kernel';

%multiref filter Values
% -expert Values for multiref 
multi_filter_values=2;

%% Apply filefilter

%do apply filefilter
% -switch for filtering the pickList by given filefilter
% -remove of carbon by polygons and bad images
do_applyFileFilter=0;

%filefilter path
% -path to the filefilter mat file (check: tom_av2_particlepicker ) 
filefilter_path='../pick/filefilter_poly.mat';


%% Extract 4 Sorting

%execute particle extraction and alignment
% -use 0 2 switch off 
do_extract_stack4sort=1;

%sort size
% -final size of extracted alg stack
sort_sz=[64 64];

%sort pixelsize (pixelsize of the template check: sort_alg_templ_path)
% -pixelsize of sort stack
sort_pixelsize=8.8;

%sort radius
% -radius of the particles on org micrograph
sort_radius=round((sort_sz(1).*(sort_pixelsize./pixelsize_micrograph))./2 );

%sort norm
% -norming for sort stack
sort_norm=['local_' num2str(round(sort_radius./2)) '&mean0+1std'];

%sort invert
% -invert flag for sorting
sort_invert=invert;

%sort invert
% -expert max size of sort stack
sort_max_sz=250000;

%sort alg tmpl path
% -path of the template for 2 2d alignment must have the same size as sort stack
sort_alg_templ_path='/fs/pool/pool-bmsan-apps/app_data/tom/autopick/26S/tmpl_cryo2d_8.8A.em';

%sort alg filter 
%values for filtering 2d alg low high smooth '' means no filter
sort_alg_filter=''; % e. g. [3 20 5]

%sort mask path
% path of the mask applied 2 sort stack
sort_mask4cl='/fs/pool/pool-bmsan-apps/app_data/tom/autopick/26S/mask_rect_62_35_8.8A.em';



%% Paralell Settings

%paralell flavour
% -use local 26Smsp, or none (switch off) 
paralell_flavour='local';

%number of processors
% -number of proc 2 use
num_proc=5;


%% Code 

if (strcmp(which('call_pick'),'/fs/pool/pool-bmsan-apps/tom_dev/av2/call_pick.m')==1)
    disp('error: you are useing the script from apps folder');
    error('copy the script 2 your working dir and execute from there!');
end;

if (strcmp(paralell_flavour,'local') || strcmp(paralell_flavour,'26Smsp'))
    try
        matlabpool close;
    catch Me
    end;
    matlabpool(paralell_flavour,num_proc(1));
end;

warning off;  mkdir(outputfolder);  warning on;
scrt_path=which('call_pick');
[a b]=unix(['cp ' scrt_path ' ' outputfolder '/call_pick_backup.m']);

if (a > 0)
    disp(b);
end;

%call picker 
if (do_pick > 0) 
   
    %transfer values
    
    %4input data
    train_st.micrograph.pixelsize=pixelsize_micrograph;
    train_st.micrograph.invert=invert;
    
    %4Matching
    train_st.template.dir=path_match_lib;
    train_st.template.pixelsize=pixelsize_match_lib;
    d=dir([path_match_lib '/*.em']);
    tmp=tom_emread([path_match_lib '/' d(1).name]);
    train_st.template.size=size(tmp.Value);
    
    
    %4Peak extraction 
    train_st.peaks.f_xcf_opt=xcf;
    train_st.peaks.f_soc_opt=soc;
    train_st.peaks.f_psr_opt=psr;
    train_st.peaks.num_of_peaks=max_num_of_particles;
    train_st.peaks.cc_thresh=cc_threshold;
    train_st.template.zero_rad=zero_rad;
    
    %4 Multiref alg
    load(mulit_path_lib);
    train_st.mult_cl=mult_cl;
    train_st.mult_cl.do_multiref=do_multiref;
    if (strcmp(multi_use_refs,'only_good'))
        train_st.mult_cl.num_cl_bad=0;
        train_st.mult_cl.cents=train_st.mult_cl.cents(:,:,1:train_st.mult_cl.num_cl_good);
    end;
    train_st.mult_cl.cc_opt=multi_cc_threshold;
    train_st.mult_cl.filter_type=multi_filter_type;
    train_st.mult_cl.filter_values=multi_filter_values;
    train_st.mult_cl.alg_max_sh=alg_max_sh;
    
    %save train_st
    warning off; mkdir([outputfolder '/pickLists']); warning on;
    save([outputfolder '/train_st.mat'],'train_st');
    
    if (iscell(img_path))
        base_p=fileparts(img_path{1});
    else
        base_p=fileparts(img_path);
    end;
    
    if (do_pick==2)
        disp('Picking Restarted: ');
        if (iscell(img_path))
            for i_t=1:length(img_path)
                [a b c]=fileparts(img_path{i_t});
                list_all(i_t).name=[b c];
            end;
        else
            list_all=dir(img_path);
        end;
        for i=1:length(list_all)
            num_all(i)=tom_get_max_numeric({list_all(i).name});
            abs_list_all{i}=[base_p '/' list_all(i).name];
        end;
        list_picked=dir([outputfolder '/pickLists/pick_*.em.mat']);
        for i=1:length(list_picked)
            num_picked(i)=tom_get_max_numeric({list_all(i).name});
        end;
        tmp_idx=find(ismember(num_all,num_picked)==0);
        disp([ num2str(length(num_picked)) ' of ' num2str(length(list_all)) ' are already picked ' ] );
        disp(['New img list has a length of: ' num2str(length(tmp_idx)) ]);
        tom_av2_picker([outputfolder '/train_st.mat'],abs_list_all(tmp_idx),[outputfolder '/pickLists']);
    else
        tom_av2_picker([outputfolder '/train_st.mat'],img_path,[outputfolder '/pickLists']);
    end;
    tom_av2_cat_pickLists([outputfolder '/pickLists/pick_*.mat'],[outputfolder '/pick_raw.mat']);
end;

%call pickList filter
if (do_applyFileFilter==1)
    unix(['cp ' outputfolder '/pick_raw.mat ' outputfolder '/pickb4filt.mat']);
    tom_av2_app_filefilter2picklist([outputfolder '/pick_raw.mat'],filefilter_path,[outputfolder '/pick_raw.mat']);
end;


%call extract 
if (do_extract_stack4sort==1)
    warning off; mkdir([outputfolder '/sort']); warning on;
    load([outputfolder '/pick_raw.mat']);
    org_alg=align2d;
    packages=tom_calc_packages(ceil(size(align2d,2)./sort_max_sz),size(align2d,2));
    
    for sp_cou=1:size(packages,1)
        align2d=org_alg(1,packages(sp_cou,1):packages(sp_cou,2));
        save([outputfolder '/sort/pick_auto_' num2str(sp_cou) '.mat'],'align2d');
        stack=tom_av2_xmipp_picklist2stack([outputfolder '/sort/pick_auto_' num2str(sp_cou) '.mat'],'','',sort_radius,sort_norm,sort_invert,sort_sz);
        template=tom_emread(sort_alg_templ_path);
        [stack all_tmpl cc_out cc_sum tmpl_nr]=tom_os3_alignStack2(stack,template.Value,'default',[0 0],'mean0+1std',sort_alg_filter,[1 5],'default');
        mask=tom_emread(sort_mask4cl);
        for i=1:size(stack,3)
            stack(:,:,i)=stack(:,:,i).*mask.Value;
        end;
        tom_emwrite([outputfolder '/sort/pick_auto_' num2str(sp_cou) '.em'],stack);
        tmpl_nr_unique=unique(tmpl_nr);
        if (size(template.Value,3) > 1)
            align2d_all_ref=align2d;
            for ii=1:length(tmpl_nr_unique)
                sub_fold=[outputfolder '/sort/ref' num2str(tmpl_nr_unique(ii)) ];
                warning off; mkdir(sub_fold); warning on;
                idx_cl=find(tmpl_nr==tmpl_nr_unique(ii));
                stack_tmp=stack(:,:,idx_cl);
                align2d=align2d_all_ref(1,idx_cl);
                tom_emwrite([sub_fold '/pick_auto_' num2str(sp_cou) '.em'],stack_tmp);
                save([sub_fold '/pick_auto_' num2str(sp_cou) '.mat'],'align2d');
            end;
        end;
    end;
    
    if (size(template.Value,3) > 1)
        cat_sub_stacks(outputfolder,size(template.Value,3),sort_max_sz);
    end;
end;


function cat_sub_stacks(outputfolder,num_of_tmpl,sort_max_sz)


for ii=1:num_of_tmpl
    fold=[outputfolder '/sort/ref' num2str(ii) ];
    dd=dir([fold '/*.em']);
    for iii=1:length(dd)
        h=tom_reademheader([fold '/' dd(iii).name]);
        sz_v(iii)=h.Header.Size(3);
    end;
    if (length(sz_v)>1)
        stack=[];
        all_alg=[];
        for isz=1:length(sz_v)
            load(strrep([fold '/' dd(isz).name],'.em','.mat'));
            tmp=tom_emreadc([fold '/' dd(isz).name]);
            stack=cat(3,stack,tmp.Value);
            all_alg=cat(2,all_alg,align2d);
        end;
        call=['rm ' fold '/*.em'];
        unix(call);
        call=['rm ' fold '/*.mat'];
        unix(call);
        packages=tom_calc_packages(ceil(sum(sz_v)./sort_max_sz),sum(sz_v));
        for ip=1:size(packages,1)
            align2d=all_alg(1,packages(ip,1):packages(ip,2));
            tmp=stack(:,:,packages(ip,1):packages(ip,2));
            tom_emwrite([fold '/pick_auto_' num2str(ip) '.em'],tmp);
            save([fold '/pick_auto_' num2str(ip) '.mat'],'align2d');
        end;
    end;
    
end;












