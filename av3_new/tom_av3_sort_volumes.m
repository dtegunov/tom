function new_idx=tom_av3_sort_volumes(p_vols,mask_path,outpath,filter,lead_z,first_vol,start_nr_out,logfile_name)
%TOM_AV3_SORT_VOLUMES sort volumes 
%
%   new_idx=tom_av3_sort_volumes(p_vols,mask_path,outpath)
%
%PARAMETERS
%
%  INPUT
%   p_vols                 path of the volumes
%   mask_path              mask for sorting
%   outpath                new path for sorted volumes
%   filter                 (0) filter kernel size     
%   lead_z                 (2) leading zeros flag 
%                              use 1 for leading zeros filenumber e.g mod_0001.em
%                              use 2 for leading zeros pre-fix e.g. 001_mod_1.em   
%   first_vol              (1 vol in the list) name of the first volume 
%   start_nr_out           (0) starting number for sorted volumes 
%   logfile_name           (outpath/_sort.log) name of the logfile
%   
%                          
%  OUTPUT
%   new_idx               sorted idx
%
%
%EXAMPLE
%
%  %example 4 filename with leading zeros
%  tom_av3_sort_volumes('out_pmbe4/models/model_iter_27*.em','out_pmbe4/mask.em','sorted/mod_',1);
%
%  %example 4 prefix numbering and defined start volume
%  tom_av3_sort_volumes('test/sp*.spi','','test_sorted/',0,2,'sp_8.spi');
%
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by FB 01/24/10
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


if (nargin < 4)
    filter=0;
end;

if (nargin < 5)
    lead_z=2;
end;

if (nargin < 6)
    first_vol='';
end;

if (nargin < 7)
    start_nr_out=1;
end;

if (nargin < 8)
    log_file_name='sort.log';
end;



d=dir(p_vols);
if (isempty(d))
    disp(['check path of input volumes']);
    error(['no volumes 4 ' p_vols]);
end;

[a b c]=fileparts(p_vols);



%build list
for i=1:length(d)
    list.path{i}=[a '/' d(i).name];
    list.used(i)=0;
end;
if (isempty(mask_path)==0)
    if (tom_isemfile(mask_path))
        mask=tom_emread(mask_path);
    else
        mask=tom_spiderread(mask_path);
    end;
    mask=mask.Value;
else
    tmp=tom_emread(list.path{1});
    mask=ones(size(tmp.Value));
end;

if (isempty(first_vol)==1)
    first=1;
else
    first=-1;
    for i=1:length(list.path)
        if (isempty(strfind(list.path{i},first_vol))==0)
            first=i;
        end;
    end;
    if (first==-1)
        error(['cannot find ' first_vol ]);
    end;
end;
list.used(first)=1;
new_idx(1)=first;


[aa bb cc]=fileparts(outpath);
% warning off;
% [au bu]=unix(['rm -r ' aa]);
% mkdir(aa);
warning on;

fid=fopen([outpath log_file_name],'wt');
disp(['using ' list.path{first} ' as first volume']);
fprintf(fid,'%s\n',['using ' list.path{first} ' as first volume']);

for i=1:length(d)-1
    
    if (tom_isemfile(list.path{first}))
        first_vol=tom_emread(list.path{first});
    else
        first_vol=tom_spiderread(list.path{first});
    end;
    
    first_vol=first_vol.Value;
    idx_unused=find(list.used(:)==0);
    cc=zeros(length(list.path),1);
    for ii=1:length(idx_unused)
        if (tom_isemfile(list.path{first}))
            part=tom_emread(list.path{idx_unused(ii)});
        else
            part=tom_spiderread(list.path{idx_unused(ii)}); 
        end;
        
        part=part.Value;
        tmp=tom_corr(first_vol,part,'norm',mask);
        [pos cc(idx_unused(ii))]=tom_peak(tmp);
     end;
    [val pos]=max(cc);
    new_idx(i+1)=pos;
    list.used(pos)=1;
    first=pos;
    disp(['Reading ' list.path{first} ]);
end;



for i=1:length(new_idx)
    if (tom_isemfile(list.path{new_idx(i)}))
        im=tom_emread(list.path{new_idx(i)});
        em_flag=1;
    else
        im=tom_spiderread(list.path{new_idx(i)});
        em_flag=0;
    end;
    
    if (filter >0)
        im.Value=tom_filter(im.Value,2);
    end;
    
    if (lead_z==0)
        if (em_flag==1)
            tom_emwrite([outpath num2str(i+start_nr_out-1) '.em'],im.Value);
            out_name_tmp=[outpath num2str(i+start_nr_out-1) '.em'];
        else
            tom_spiderwrite([outpath num2str(i+start_nr_out-1) '.spi'],im.Value);
            out_name_tmp=[outpath num2str(i+start_nr_out-1) '.spi'];
        end;
    end;
    if (lead_z==1)
         if (em_flag==1)
            tom_emwrite([outpath get_zeros(i-1+start_nr_out) '.em'],im.Value);
            out_name_tmp=[outpath get_zeros(i-1+start_nr_out) '.em'];
         else
            tom_spiderwrite([outpath get_zeros(i-1+start_nr_out) '.spi'],im.Value);
            out_name_tmp=[outpath get_zeros(i-1+start_nr_out) '.spi'];
         end;
    end;
    if (lead_z==2)
         [aa bb cc]=fileparts(outpath);
         [aa_out bb_out cc_out]=fileparts(list.path{new_idx(i)});
         if (em_flag==1)
            tom_emwrite([aa '/' get_zeros(i-1+start_nr_out) bb_out cc_out],im.Value);
         else
            tom_spiderwrite([aa '/' get_zeros(i-1+start_nr_out) bb_out cc_out],im.Value);
         end;
         out_name_tmp=[aa '/' get_zeros(i-1+start_nr_out) bb_out cc_out];
    end;
    
    fprintf(fid,'%s\n',[list.path{new_idx(i)} ' ==> ' out_name_tmp]);
    disp([list.path{new_idx(i)} ' ==> ' out_name_tmp]);
end;
fclose(fid);

function out_str=get_zeros(num)

if (num<10)
    my_zer='0000';
end;

if (num>=10 && num<100)
    my_zer='000';
end;

if (num>=100 && num<1000)
    my_zer='00';
end;

if (num>=1000 && num<10000)
    my_zer='00';
end;

if (num>=10000 && num < 100000)
    my_zer='0';
end;

if (num>=100000)
    my_zer='';
end;

out_str=[my_zer num2str(num)]; 


