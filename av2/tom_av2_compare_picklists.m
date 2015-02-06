function [total eval merged_ref_test al_over al_fp al_fn merged_ov_fp_fn]=tom_av2_compare_picklists(ref,test,mode,radius,find_what,replace_with,hd_output,verbose)
%TOM_AV2_COMPARE_PICKLISTS compares 2 picklists
%
%   eval=tom_av2_compare_picklists(hand,test)
%
%PARAMETERS
%
%  INPUT
%   ref                 reference picklist
%   test                test picklist will be compared against ref
%   mode                ('gaps') or 'no_gaps'
%                        gaps:      only images which appear in both lists are used   
%                        no_gaps:   all images in test and ref are used.
%                                   (e.g. if there is no correspondig image in test all parts of test are false psitive) 
%   radius               (ref(1,1).radius) search radius
%   find_what            rpl string to adapt test path 2 ref path or use DObasename for matching filenames only 
%   replace_with         rpl string to adapt test path 2 ref path  
%   hd_output            (1) write combined picklists 2 hardisk use 0 toswitch off
%   verbose              vebose flag 
%
%  OUTPUT
%  total               containing the total
%  eval                eval struct containing the result per imgae
%  merged_ref_test     merged picklist for tom_av2_particlepickergui 
%  al_over             picklist with overlap 
%  al_fp               picklist containing false positive
%  al_fn               picklist containing false negative 
%  merged_ov_fp_fn     merged picklist for tom_av2_particlepickergui 
%
%EXAMPLE
%    
%    [total eval]=tom_av2_compare_picklists('pl_fb.mat','align-26S_1003.em.mat','gaps',80,'/fs/sun16/lv01/pool/pool-nickell3/fb/4Tueb/271206/high_corr/','/fs/sun16/lv01/pool/pool-nickell3/fb/4Hrabe/inter_picker_var/high/',1);
%    
%    [total eval align2d]=tom_av2_compare_picklists('pl_fb.mat','align-26S_1003.em.mat','gaps',80,'DObasename');
%    
%     %without hd output 
%    [total eval]=tom_av2_compare_picklists(hand,'out/align-p47f11_3860.em.mat','gaps',60,'DObasename','',0); 
%
%
%REFERENCES
%
%SEE ALSO
%   ...
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

%load listst 

if (nargin<3)
    mode='gaps';
end;

if (nargin<4)
    radius=ref(1,1).radius;
end;

if (nargin<5)
    find_what='';
end;

if (nargin<6)
    replace_with='';
end;

if (nargin<7)
    hd_output=1;
end

if (nargin<8)
    verbose=1;
end

if (verbose==1)
    disp_param(ref,test,mode,radius,find_what,replace_with,verbose);
end;

if (isstruct(ref)==0)
    load(ref);
    ref=align2d;
end;




if (isstruct(test)==0)
    if exist(test,'file')
        load(test);
        test=align2d;
        clear('align2d');
    else
        dd=dir(test);
        bp=fileparts(test);
        disp(['reading ' num2str(length(dd)) ' files']);
        if (isempty(bp))
            load([dd(1).name]);
        else
            load([bp '/' dd(1).name]);
        end;
        
        test=align2d;
        for i=2:length(dd)
            if (isempty(bp))
                load([dd(i).name]);
            else
                load([bp '/' dd(i).name]);
            end;
            test=cat(2,test,align2d);
            if (mod(i,100)==0)
                disp([num2str(i) ' lists loaded']);
            end;
        end;
        clear('align2d');
    end;
end;

for i=1:size(ref,2)
    if (strcmp(find_what,'DObasename'))
        [all_bases{i} b c]=fileparts(ref(1,i).filename);
        ref(1,i).filename=[b c];
    end;
end;

u_base_p=unique(all_bases);

if (length(u_base_p)>1)
    error('multible basepath found');
end;


%get unique images
for i=1:size(ref,2)
    all_ref_names{i}=ref(1,i).filename;
end;
u_ref_names=unique(all_ref_names);


for i=1:size(test,2)
    if (strcmp(find_what,'DObasename'))
        [a b c]=fileparts(test(1,i).filename);
        test(1,i).filename=[b c];
    else
        test(1,i).filename=strrep(test(1,i).filename,find_what,replace_with);
    end;
end;

for i=1:size(test,2)
    all_test_names{i}=test(1,i).filename;
end;
u_test_names=unique(all_test_names);

if (strcmp(mode,'gaps'))
    both_names=intersect(u_ref_names,u_test_names);
    if (isempty(both_names))
        disp('no corresponding images found ...try find_what replace_with parameter');
        disp(['ref img 1: ' u_ref_names{1}]);
        disp(['test img 1: ' u_test_names{1}]);
        return;
    else
        if (verbose==1)
            disp([num2str(length(both_names)) ' images in both lists found']);
        end;
    end;
else
    both_names=unique(union(u_ref_names,u_ref_names));
end;

%alloc memory
names={};
img_idx=zeros(length(both_names),1);
p_nr_test=zeros(length(both_names),1);
p_nr_ref=zeros(length(both_names),1);
fp=zeros(length(both_names),1);
fph=zeros(length(both_names),1);
fn=zeros(length(both_names),1);
overlap=zeros(length(both_names),1);
noFound=zeros(length(both_names),1);
falsePositives=zeros(length(both_names),1);

parfor i=1:length(both_names)
%for i=1:length(both_names)
    %find hits in both
    names{i}=both_names{i};
    [aa bb c]=fileparts(both_names{i});
    [a b]=strtok(bb,'_');
    img_idx(i)=str2double(strrep(b,'_',''));
    tmp_ref=ref(1,find(ismember(all_ref_names,both_names{i})));
    tmp_test=test(1,find(ismember(all_test_names,both_names{i})));
    p_nr_test(i)=size(tmp_test,2);
    p_nr_ref(i)=size(tmp_ref,2);
    [noFound(i) falsePositives(i) overlap(i) fp(i) fph(i) fn(i) al_over_tmp{i} al_fp_tmp{i} al_fn_tmp{i}] = comparePicks(tmp_test,tmp_ref,radius-1);
    if (verbose==1)
        disp([bb c ' num ref: ' num2str(size(tmp_ref,2))  ' num test: ' num2str(size(tmp_test,2)) ' overlap : ' num2str(overlap(i)) ' false positive: ' num2str(fp(i)) ' false negative: ' num2str(fn(i)) ]);
    end;
end;
    
eval.names=names;
eval.img_idx=img_idx;
eval.p_nr_test=p_nr_test;
eval.p_nr_ref=p_nr_ref;
eval.fp=fp;
eval.fph=fph;
eval.fn=fn;
eval.overlap=overlap;

total.p_nr_test=sum(p_nr_test);
total.p_nr_ref=sum(p_nr_ref);
total.fp=sum(falsePositives)./size(test,2);
total.fph=sum(falsePositives)./total.p_nr_ref;
total.fn=1-(sum(noFound)./total.p_nr_ref);
total.overlap=sum(noFound)./total.p_nr_ref;


if (verbose==1)
    disp(' ');
    disp('Total: ');
    disp(['num ref: ' num2str(total.p_nr_ref)  ' num test: ' num2str(total.p_nr_test) ' overlap : ' num2str(total.overlap) ' false positive: ' num2str(total.fp) ' false negative: ' num2str(total.fn) ]);
    disp(['num overlap: ' num2str(sum(noFound))  ' num false positive: ' num2str(sum(falsePositives)) ' num false negative ' num2str(total.p_nr_ref-sum(noFound)) ]);   
end;

if (nargout>2 || hd_output==1)
    
    [merged_ref_test merged_ov_fp_fn al_over al_fp al_fn ]=transform_and_merge(ref,test,al_over_tmp,al_fp_tmp,al_fn_tmp,u_base_p{1},ref(1,1).radius,find_what,hd_output);
    
end;

function [merged_ref_test merged_ov_fp_fn al_over al_fp al_fn ]=transform_and_merge(ref,test,al_over_tmp,al_fp_tmp,al_fn_tmp,base_p,rad,find_what,hd_output)


if strcmp(find_what,'DObasename')
    base_flag=1;
else
    base_flag=0;
end;

%transform
al_over=al_over_tmp{1};
al_fp=al_fp_tmp{1};
al_fn=al_fn_tmp{1};

for i=2:length(al_over_tmp)
    al_over=cat(2,al_over,al_over_tmp{i});
    al_fp=cat(2,al_fp,al_fp_tmp{i});
    al_fn=cat(2,al_fn,al_fn_tmp{i});
end;

%clean and extend path
for i=1:size(al_over,2)
    if (base_flag==1)
        tmp(1,i).filename=[base_p '/' al_over(1,i).filename];
    else
        tmp(1,i).filename=[al_over(1,i).filename];
    end;
    tmp(1,i).radius=rad;
    tmp(1,i).position.x=al_over(1,i).position.x;
    tmp(1,i).position.y=al_over(1,i).position.y;
    tmp(1,i).ccc=al_over(1,i).ccc;
    tmp(1,i).class='overlap';
    tmp(1,i).color=[0 1 0];
end;
if (exist('tmp','var'))
    al_over=tmp;
else
    al_over=[];
end;

clear('tmp');

for i=1:size(al_fp,2)
    if (base_flag==1)
        tmp(1,i).filename=[base_p '/' al_fp(1,i).filename];
    else
        tmp(1,i).filename=[al_fp(1,i).filename];
    end;
    tmp(1,i).radius=rad;
    tmp(1,i).position.x=al_fp(1,i).position.x;
    tmp(1,i).position.y=al_fp(1,i).position.y;
    tmp(1,i).ccc=al_fp(1,i).ccc;
    tmp(1,i).class='falsePositive';
    tmp(1,i).color=[1 0 0];
end;
if (exist('tmp','var'))
    al_fp=tmp;
else
    al_fp=[];
end;

clear('tmp');

for i=1:size(al_fn,2)
    if (base_flag==1)
        tmp(1,i).filename=[base_p '/' al_fn(1,i).filename];
    else
        tmp(1,i).filename=[al_fn(1,i).filename];
    end;
    tmp(1,i).radius=rad;
    tmp(1,i).position.x=al_fn(1,i).position.x;
    tmp(1,i).position.y=al_fn(1,i).position.y;
    tmp(1,i).ccc=al_fn(1,i).ccc;
    tmp(1,i).class='falseNegative';
    tmp(1,i).color=[0 0 1];
end;
if (exist('tmp','var'))
    al_fn=tmp;
else
    al_fn=[];
end;

if (isempty(al_over) &&  isempty(al_fp) && isempty(al_fn))
    merged_ov_fp_fn=[];
end;

if (isempty(al_over) &&  isempty(al_fp) && isempty(al_fn)==0 )
    merged_ov_fp_fn=al_fn;
end;

if (isempty(al_over) &&  isempty(al_fp)==0 && isempty(al_fn) )
    merged_ov_fp_fn=al_fp;
end;

if (isempty(al_over) &&  isempty(al_fp)==0 && isempty(al_fn)==0 )
    merged_ov_fp_fn=cat(2,al_fp,al_fn);
end;

if (isempty(al_over)==0 &&  isempty(al_fp) && isempty(al_fn) )
    merged_ov_fp_fn=al_over;
end;

if (isempty(al_over)==0 &&  isempty(al_fp) && isempty(al_fn)==0 )
    merged_ov_fp_fn=cat(2,al_over,al_fn);
end;

if (isempty(al_over)==0 &&  isempty(al_fp)==0 && isempty(al_fn))
    merged_ov_fp_fn=cat(2,al_over,al_fp);
end;

if (isempty(al_over)==0 &&  isempty(al_fp)==0 && isempty(al_fn)==0)
    merged_ov_fp_fn=cat(2,al_over,al_fp,al_fn);
end;



zz=1;
rad=ref(1,1).radius;
for i=1:size(ref,2)
    merged_ref_test(1,zz).filename=ref(1,i).filename;
    merged_ref_test(1,zz).radius=rad;
    merged_ref_test(1,zz).position.x=ref(1,i).position.x;
    merged_ref_test(1,zz).position.y=ref(1,i).position.y;
    merged_ref_test(1,zz).ccc=ref(1,i).ccc;
    merged_ref_test(1,zz).class='ref';
    merged_ref_test(1,zz).color=[0 1 0];
    zz=zz+1;
end;
for i=1:size(test,2)
    merged_ref_test(1,zz).filename=test(1,i).filename;
    merged_ref_test(1,zz).radius=rad;
    merged_ref_test(1,zz).position.x=test(1,i).position.x;
    merged_ref_test(1,zz).position.y=test(1,i).position.y;
    merged_ref_test(1,zz).ccc=test(1,i).ccc;
    merged_ref_test(1,zz).class='test';
    merged_ref_test(1,zz).color=[1 0 0];
    zz=zz+1;
end;

if (hd_output==1)
    
    try
        warning off;
        mkdir('comb_results')
        warning on;
        align2d=merged_ref_test;
        save('comb_results/merged_ref_test.mat','align2d');
        align2d=merged_ov_fp_fn;
        save('comb_results/merged_ov_fp_fn.mat','align2d');
        align2d=al_fn;
        save('comb_results/al_fn.mat','align2d');
        align2d=al_fp;
        save('comb_results/al_fp.mat','align2d');
        align2d=al_over;
        save('comb_results/al_over.mat','align2d');
    catch Me
        disp('error writing results 2 disk');
        disp(Me.message);
    end;
    
    
end;


function [noFound falsePositives pF pFP pFPh pFN al_ov al_fp al_fn] = comparePicks(al1,al2,radius)

%al2 == hand ground trouth
%al1 == auto test set

noFound = 0;
falsePositives = 0;
currentFile = al1(1).filename;

al2_org=al2;

for a1_iter=1:size(al1,2)
    
    if(~strcmp(currentFile,al1(a1_iter).filename))
        currentFile = al1(a1_iter).filename;
    end;
    
    [found idx_al2]= findInAl(al1(1,a1_iter),al2,currentFile,radius);
    
    if(found)
        noFound = noFound +1;
        al_ov(1,noFound)=al1(1,a1_iter);
        al2_found(noFound)=idx_al2;
        al2(1,idx_al2).position.x=-2*radius-1; % remove from search
        al2(1,idx_al2).position.y=-2*radius-1; % remove from search
    else
        falsePositives = falsePositives +1;
        al_fp(1,falsePositives)=al1(1,a1_iter);
    end;
    
end;
al2=al2_org;

%get al false negative
if (noFound==0)
    al2_found=-1;
end;
idx=setdiff(1:size(al2,2),al2_found);
al_fn=al2(1,idx);

if (exist('al_ov','var')==0)
    al_ov=[];
end;
if (exist('al_fp','var')==0)
    al_fp=[];
end;

if (exist('al_fn','var')==0)
    al_fn=[];
end;


pF = noFound / length(al2);
pFP = falsePositives / length(al1);
pFPh = falsePositives / length(al2);
pFN = 1 - pF;

function [found found_idx]= findInAl(pick,al2,currentFile,radius)

found = false;

a2_iter = 0;
pos1 = [pick.position.x pick.position.y];
found_idx=[];

all_pos_al2=zeros(size(al2,2),2);
for i=1:size(al2,2)
    all_pos_al2(i,1)=al2(1,i).position.x;
    all_pos_al2(i,2)=al2(1,i).position.y;
end;
 
[found_idx, f_coord, dist] = tom_nearestpoint(pos1,all_pos_al2);
if (dist <= radius)
    found=1;
else
    found_idx=[];
end;





function disp_param(ref,test,mode,radius,find_what,replace_with,verbose)

disp(' ');
disp('================>');
if (isstruct(ref)==0)
    disp(['ref: ' ref]);
else
    disp(['ref: struct in mem']);
end;
    
if (isstruct(test)==0)
    disp(['test: ' test]);
else
    disp(['test: struct in mem']);
end;

disp(['mode: ' mode]);
disp(['radius: ' num2str(radius)]);
disp(['find_what: ' find_what]);
disp(['replace_with: ' replace_with]);
disp(['verbose: ' num2str(verbose)]);
disp('<================');
disp(' ');






























