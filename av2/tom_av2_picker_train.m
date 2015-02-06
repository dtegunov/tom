function train_st=tom_av2_picker_train(model3d,picklist,outputfolder,options)
%tom_av2_picker_train trains a 2d picker
%
%   train_st=tom_av2_picker_train(model3d,picklist,options)
%
%PARAMETERS
%
%  INPUT
%   model3d             3d model (pixelsize must be the same as in the picklist)
%   picklist            picklist in tom format 
%   outputfolder        folder for storing projections 
%   options             contains all options (... for experts 2 guide the training)
%
%
%  OUTPUT
% 
%   train_st            training structure
%
%EXAMPLE
%    
% matlabpool open local 8;
% train_st=tom_av2_picker_train('model_rot.em','train.mat','output3');
%
%REFERENCES
%
%SEE ALSO
%   tom_av2_picker
%
%   created by SN/FB 01/24/06
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




disp(['Reading ' picklist]);
load(picklist);
hand=align2d; clear('align2d');
[all test al_test train al_train num_p num_img]=extract_pickList(hand,36);
disp(['picklist contains: ' num2str(num_p) ' particles  on ' num2str(num_img) ' images.']);

disp(['Reading ' model3d]);
mod=tom_emread(model3d);
mod=mod.Value;
disp(' ');

if (nargin<4)
    options=gen_default_options(outputfolder,mod);
end;



disp(['Projection model (bin) ' num2str(options.template.bin)]);
proj=gen_proj(tom_bin(mod,options.template.bin),options.template.proj_ang.phi,options.template.proj_ang.psi,options.template.proj_ang.theta);

disp(['Compressing projections']);
comp_proj=tom_os3_classify_templStack(proj,'',options.template.num_cl,1,1,'classify',0,0);
clear('proj');

disp(['Writing projections 2 ' options.template.dir]);
write_proj(comp_proj,options.template.dir);


disp(['Starting filter optimization']);
warning off;
    mkdir([outputfolder '/filter_train']);
    mkdir([outputfolder '/filter_test']);
warning on;


for ii=1:length(options.peaks.f_xcf)
    parfor i=1:length(train)
    %for i=1:length(train) 
        [num(i) fn(i)]=opt_match(train{i},hand,options);
    end;
    all_f_xcf(ii)=options.peaks.f_xcf(ii);
    all_num_mean(ii)=round(mean(num));
    all_fn_mean(ii)=mean(fn);
    disp(['iter: ' num2str(ii) ' f_xcf: ' num2str(all_f_xcf(ii)) ' mean_num: '  num2str( all_num_mean(ii)) ]);
    idx_val=find(fn>options.peaks.max_part_lost);
    if (isempty(idx_val)==0)
        disp(['Warning: ' num2str(length(idx_val)) ' of '  num2str(length(num)) ' have the specified part lost' ]);
        disp(['mean fn is: ' num2str(mean(all_fn_mean(ii))) ]);
    end;
 
end;
[min_val min_pos]=min(all_num_mean);
min_pos=min(min_pos);


options.peaks.f_xcf_opt=all_f_xcf(min_pos);
options.peaks.f_soc_opt=1-all_f_xcf(min_pos);
options.peaks.num_of_peaks=all_num_mean(min_pos);

clear('fn');
disp(['Testing with: num. of peaks: ' num2str(options.peaks.num_of_peaks) ' f_xcf: ' num2str(options.peaks.f_xcf)  '  f_soc ' num2str(options.peaks.f_soc)  ]);
parfor i=1:length(test)
     [all_aling{i} fn(i) fp(i) over(i)]=match(test{i},hand,options);
end;

disp(['Rates: fn: ' num2str(mean(fn)) ' fp ' num2str(mean(fp)) ' overlap: ' num2str(mean(over)) ]);

disp(['Generating Particle stack 4 classifier optimization']);
[all test al_test train al_train num_p num_img]=extract_pickList(hand,24);

parfor i=1:length(train)
%for i=1:length(train)
    all_aling{i}=match(train{i},hand,options);
end;

al_filter=all_aling{1};
for i=2:size(all_aling,2)
    al_filter=cat(2,al_filter,all_aling{i});
end;


rad=round((options.template.size(1).*2^options.template.bin)./2);
[total eval merge al_over al_fp al_fn]=tom_av2_compare_picklists(hand,al_filter,'gaps',rad,'DObasename','',0,0); 

stack_bad=tom_av2_xmipp_picklist2stack(al_fp,'','',rad,'gradient&mean0+1std',1,options.mult_cl.bin);
stack_good=tom_av2_xmipp_picklist2stack(al_train,'','',rad,'gradient&mean0+1std',1,options.mult_cl.bin);


disp(['Starting Classifyer training']);

options.test_good_bad_ratio=(num_p./num_img)./options.peaks.num_of_peaks;

warning off;
    mkdir([outputfolder '/classifyer_train']);
warning on;

tom_emwrite([outputfolder '/stack_good.em'],stack_good);
tom_emwrite([outputfolder '/stack_bad.em'],stack_bad);

mult_cl=tom_av2_cl_multiref_train(stack_good,stack_bad,options.mult_cl);

options.mult_cl=mult_cl;
train_st=options;

save([outputfolder '/train_st.mat'],'train_st');

warning off;
mkdir([outputfolder '/pick_test']);
warning on;

tom_av2_picker(train_st,test,[outputfolder '/pick_test']);


[total eval merge al_over al_fp al_fn]=tom_av2_compare_picklists(hand,[outputfolder '/pick_test/pick_*.mat'],'gaps',rad,'DObasename','',0,0); 

disp(' ');
disp('Total: ');
disp(['num ref: ' num2str(total.p_nr_ref)  ' num test: ' num2str(total.p_nr_test) ' overlap : ' num2str(total.overlap) ' false positive: ' num2str(total.fp) ' false negative: ' num2str(total.fn) ]);
%disp(['num overlap: ' num2str(sum(noFound))  ' num false positive: ' num2str(sum(falsePositives)) ' num false negative ' num2str(total.p_nr_ref-sum(noFound)) ]);


function [align2d fn_out fp_out over_out]=match(img_name,hand,options)


im=tom_emreadc(img_name);
im=tom_bin(im.Value,options.template.bin);
peaks=tom_os3_match_templ(im,options.template.dir,options.peaks);
[coord,align2d]=tom_os3_extract_peaks2(peaks,options.peaks.num_of_peaks,options.template.size,img_name,options.template.bin);
detect_rad=round((options.template.size(1).*2^options.template.bin)./2)-round((options.template.size(1)*2^options.template.bin)./2.*0.7);
[total eval]=tom_av2_compare_picklists(hand,align2d,'gaps',detect_rad,'DObasename','',0,0);

if (isempty(find(coord(:) > 4097))==0)
    disp(' ');
end;

fn_out=total.fn;
fp_out=total.fp;
over_out=total.overlap;

function [num_of_peaks_out fn_out]=opt_match(img_name,hand,options)


im=tom_emreadc(img_name);
im=tom_bin(im.Value,options.template.bin);
peaks=tom_os3_match_templ(im,options.template.dir,options.peaks);
max_num_of_peaks=options.peaks.max_num_of_peaks;

[coord,align2d]=tom_os3_extract_peaks2(peaks,max_num_of_peaks,options.template.size,img_name,options.template.bin);

detect_rad=round((options.template.size(1).*2^options.template.bin)./2)-round((options.template.size(1)*2^options.template.bin)./2.*0.7);

zz=0;
incre=15;
for i=1:incre:max_num_of_peaks
    al_cut=align2d(1,1:i);
    [total eval]=tom_av2_compare_picklists(hand,al_cut,'gaps',detect_rad,'DObasename','',0,0);
    zz=zz+1;
    fn(zz)=total.fn;
    inter_v(zz)=i;
end;
[val pos]=min(abs([fn-options.peaks.max_part_lost]));
pos=min(pos);
num_of_peaks_out=inter_v(pos);
[total eval]=tom_av2_compare_picklists(hand,align2d(1,1:num_of_peaks_out),'gaps',detect_rad,'DObasename','',0,0);
fn_out=total.fn;


    
function [all test align_test train align_train num_p num_img]=extract_pickList(align2d,max_sz)

for i=1:size(align2d,2)
    all_f{i}=align2d(1,i).filename;
end;
all=unique(all_f);

idx=randperm(length(all));

if length(idx)>max_sz
    idx=idx(1:max_sz);
end;

test=all(idx(1:round(length(idx).*0.333333)));
train=all(idx(round(length(idx).*0.333333)+1:end));

idx_align_test=find(ismember(all_f,test));
align_test=align2d(1,idx_align_test);
idx_align_train=find(ismember(all_f,train));
align_train=align2d(1,idx_align_train);

num_p=size(align2d,2);
num_img=(length(all));





function proj=gen_proj(model,phi,psi,theta)

zz=1;
for phi_iter=phi
    for psi_iter=psi
        for theta_iter=theta
            ang(:,zz)=tom_sum_rotation([0 psi_iter theta_iter; 270 90 phi_iter],[0 0 0; 0 0 0]);
            zz=zz+1;
        end;
    end;
end;

proj=zeros(size(model,1),size(model,2),size(ang,2));

parfor iii=1:zz-1
    proj(:,:,iii)=sum(tom_rotate(model,ang(:,iii)'),2);
    if (mod(iii,100)==0)
        fprintf('.'); drawnow;
    end;
end;
disp('done!');

function write_proj(comp_proj,outputfolder)


outfilestruct.method='singleFiles';
outfilestruct.num_of_entries=100000;
outfilestruct.path=outputfolder;
outfilestruct.folder_num_offset=-1;
outfilestruct.filename='template_';
outfilestruct.ext='.em';
outfilestruct.fileformat='em';
tom_av2_stack_clean(outfilestruct);
tom_av2_stack_write(comp_proj,outfilestruct,1);


function options=gen_default_options(output,model)


options.template.bin=3;
options.template.proj_ang.phi=[1:12:348];
options.template.proj_ang.psi=[1:10:60];
options.template.proj_ang.theta=[1:12:348];
options.template.num_cl=100;
%options.template.num_cl=20;
options.template.dir=[output '/templates'];
options.template.size=[round(size(model,1)./(2^options.template.bin)) round(size(model,2)./(2^options.template.bin))];

options.peaks.max_num_of_peaks=300;
options.peaks.num_of_peaks=118;
options.peaks.max_part_lost=0.08;
options.peaks.f_xcf=[0.10];
options.peaks.f_soc=[0.90];
options.peaks.f_psr=0;
options.peaks.f_xcf_opt=0.15;
options.peaks.f_soc_opt=0.85;
options.peaks.f_psr_opt=0;

%options.mult_cl.nr_good_train=[100 200 400 400];
%options.mult_cl.nr_bad_train=[400 800 1600 400];
options.mult_cl.nr_good_train=[100 300 ];
options.mult_cl.nr_bad_train= [500 1500];
options.mult_cl.eigs=[1 20];
options.mult_cl.bin=3;
options.mult_cl.num_cl_good=[1 3 3 3];
options.mult_cl.num_cl_bad=[1 3 5 6];
% options.mult_cl.num_cl_good=[1 3 6];
% options.mult_cl.num_cl_bad=[1 3 6];
options.mult_cl.num_al=5;
options.mult_cl.cc_search_range=0.2:0.02:0.42;
options.mult_cl.test_good_bad_ratio=0.25;
options.mult_cl.test_max_num=1000;
 


