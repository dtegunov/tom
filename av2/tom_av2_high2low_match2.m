function new_idx=tom_av2_high2low_match2(f_align_high,sel_high,f_align_low,sel_low,post_fix,max_sh,max_dist,bin)
%TOM_AV2_HIGH2LOW_MATCH2 matches a high to low picklist
%
%   tom_av2_high2low_match2(f_align_high,sel_high,f_align_low,sel_low,post_fix,max_sh,max_dist,bin)
%
%  TOM_AV2_HIGH2LOW_MATCH2 matches a high 2 low picklist and gives a
%  htl(high 2 low) file back
%  
%
%PARAMETERS
%
%  INPUT
%   f_align_high        filename of the high alignment file
%   sel_high            filename of the high sel file
%   f_align_low         filename of the low alignment file
%   sel_low            sel file name for the low images
%   post_fix           ('corr') post-fix for folder names 
%   max_sh             (500) max shift between the high an low images 
%   max_dist           (200) max dist for machting use particle radius
%   bin                (1) binning for calculation of the shift   
%   
%  OUTPUT
%   new_idx          gives the high2low index in memory ...just for testing
%                     and plotting
%                    
%                 ...the htl and index.mat file are written to dist with
%                 the filenmame [f_align_high f_align_low .htl]
%                        
%
%EXAMPLE
%     
%   new_idx=tom_av2_high2low_match2('11_corrf_high_128.mat','11_corrf_high_128.sel','11_corrf_low_128.mat','11_corrf_low_128.sel','corr',400,100,1);
%   
%
%   tom_av2_high2low_match2('18_corrf_low_128.mat','18_corrf_low_128.sel');
%   %example for no pairs ...just lows   
%
%REFERENCES
%
%SEE ALSO
%   tom_av2_check_htl_file
%
%   created by fb
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



if (nargin <3)
    disp('running fake htl ...cell file column is just copied !');
    disp('only lows exist !');
    name_out=[f_align_high '_' f_align_high '.htl'];
    call=['awk ''{ printf "%s %s\n",$1,$1 } '' ' sel_high ' > ' name_out];
    unix(call);
    disp(call);
    tom_av2_htl_conv2index(name_out);
    return;
end;



if (nargin < 5)
    post_fix='corr';
end;

if (nargin < 6)
    max_sh=500;
end;

if (nargin < 7)
    max_dist=200;
end;

if (nargin < 8)
    bin=1;
end;



low_name=['/low_' post_fix];
high_name=['/high_' post_fix];



%load and initialize !
load(f_align_high);
align_high=align2d;
load(f_align_low);
align_low=align2d;

new_idx=zeros(size(align_high,2),1);


%get unique list of filenames
for i=1:size(align_low,2)
    all_names_low{i}=align_low(1,i).filename;
end;
name_list_low=unique(all_names_low);
clear('all_names_low');

%get unique list of filenames
for i=1:size(align_high,2)
    all_names_high{i}=align_high(1,i).filename;
end;
name_list_high=unique(all_names_high);
clear('all_names_high');


look_count=1;
%build up lookuptable
five_p=round(length(name_list_low)./20);
zz_p=1;
tic;



for i=1:length(name_list_low)
    
    
    idx=find(ismember(name_list_high,strrep(name_list_low{i},low_name,high_name)));
    
    %idx=find(ismember(name_list_high,strrep(name_list_low{i},low_name,'/high_corr')));
    
    if (isempty(idx))
        continue;
    end;
    
    
%     if (isempty(post_fix))
%         headr=tom_emheader(name_list_low{i});
%         sh=str2num(char(headr.Header.Comment)');
%     else
%         headr=tom_emheader(strrep(name_list_low{i},'_corr/','/'));
%         sh=str2num(char(headr.Header.Comment)');
%     end;
    
    
%     if (isempty(sh)==0 &&  ndims(sh)==2)
%         shift=sh;
%         
    %else rmove comment
     
       pos_low=get_pos_per_img(align_low,name_list_low{i});
       pos_high=get_pos_per_img(align_high,strrep(name_list_low{i},low_name,high_name));
       shift=calc_shift(pos_high(:,1:2)',pos_low(:,1:2)',max_sh,bin);
   % end
    
    %disp(['Hack Warning: diff is '  num2str(sh-shift) ]);
    
    %apply shift 2 pos low
    pos_low=pos_low+repmat([shift 0],size(pos_low,1),1);
    
    %save in lookup table
    lookup.filename{look_count}=name_list_low{i};
    lookup.pos_and_id{look_count}=pos_low;
    lookup.shift{look_count}=shift;
    look_count=look_count+1;
    
    if (mod(i,five_p)==0)
        toc;
        disp([num2str(zz_p.*5) '% done'] );
        zz_p=zz_p+1;
        tic;
    end;
    
end;

match_count=0;
no_match_count=0;
match_dist=zeros(size(align_high,2),1);
match_dist=match_dist-1;

for i=1:size(align_high,2)
    idx=find(ismember(lookup.filename,strrep(align_high(1,i).filename,high_name,low_name)));
    if (isempty(idx)==0)
        pos_and_ind=lookup.pos_and_id{idx};
        [pointidx, pointcoords, distance] = tom_nearestpoint([align_high(1,i).position.x align_high(1,i).position.y],pos_and_ind(:,1:2));
        pos_and_ind(pointidx,1:2)=[-90000 -90000];
        lookup.pos_and_id{idx}=pos_and_ind;
        
        if(distance < max_dist)
            new_idx(i)=pos_and_ind(pointidx(1),3);
            match_count=match_count+1;
            match_dist(match_count)=distance;
        else
            new_idx(i)=0;
            no_match_count=no_match_count+1;
        end;
    end;
end;

tmp_match_dist=match_dist(1:(match_count));
clear('match_dist'); match_dist=tmp_match_dist; clear('tmp_match_dist');


disp(['Particles in low alignment File: ' num2str(size(align_low,2)) ] );
disp(['Particles in high alignment File: ' num2str(size(align_high,2)) ] );
disp(' ');
disp('Match dist stat:');
tom_dev(match_dist);
disp(['Particles in Range: ' num2str(match_count) '  Particles out of Range:  ' num2str(no_match_count) ] );
disp(' ');


filename_out=[f_align_high '_' f_align_low '.htl'];

fid=fopen(filename_out,'w');

sel_h=importdata(sel_high);
sel_l=importdata(sel_low);

for i=1:size(align_high,2)
    if (new_idx(i)>0)
        name_high=sel_h.textdata{i};
        name_low=sel_l.textdata{new_idx(i)};
        fprintf(fid,[name_high ' ,'  name_low ' \n']);
    end;
    
end;

fclose(fid);

save([f_align_high '_' f_align_low '.mat'],'new_idx');


[isunique lines]=tom_av2_xmipp_check_unique([f_align_high '_' f_align_low '.htl']);

disp('adding 2 columns with part idx');
tom_av2_htl_conv2index([f_align_high '_' f_align_low '.htl']);

if (isunique==0)
    error('dataset corrupted check input!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
end;

try
    unix(['chmod ugo+rwx ' f_align_high '_' f_align_low '.*' ]);
catch
    disp('Cannot chmod !');
end;


function pos=get_pos_per_img(align,img_name)


pos=zeros(40000,3);
zz=1;
for i=1:size(align,2)
    if (strcmp(align(1,i).filename,img_name))
        pos(zz,:)=[align(1,i).position.x align(1,i).position.y i];
        zz=zz+1;
    end;
end;

if (zz==1)
    pos=-1;
    return;
end;

tmp=pos(1:zz-1,:);
clear('pos');

pos=tmp;


function sh=calc_shift(pos_high,pos_low,max_sh,bin)

pos_high=round(pos_high./2^bin);
pos_low=round(pos_low./2^bin);

max_sh=max_sh./2^bin;

max_sz=max([max(max(pos_low)) max(max(pos_high))]);
mid=floor(max_sz./2)+1;
mask=tom_spheremask(ones(max_sz,max_sz),max_sh);

tmp_img_low=zeros(max_sz,max_sz);
tmp_img_high=zeros(max_sz,max_sz);

for i=1:size(pos_high,2)
    tmp_img_high(pos_high(1,i),pos_high(2,i))=1;
end;

for i=1:size(pos_low,2)
    if (pos_low(1,i) < 1 || pos_low(2,i) < 1)
        continue;
    end;
    tmp_img_low(pos_low(1,i),pos_low(2,i))=1;
end;

cc=tom_corr(tmp_img_low,tmp_img_high,'norm');

[a b]=tom_peak(cc.*mask);

a=a(1,:);

sh=a-mid;

sh=sh.*2^bin;

