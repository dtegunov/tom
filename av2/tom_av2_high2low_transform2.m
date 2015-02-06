function tom_av2_high2low_transform2(input_sel,input_doc,input_htl,output_sel,f_output_doc,dir_flag,all_flag)
%tom_av2_high2low_transform2(input_sel,input_doc,input_htl,output_sel,f_output_doc,dir_flag,all_flag)
%
%   tom_av2_high2low_transform2(input_sel,input_htl,output_sel,flag,match_flag)
%
%  TOM_AV2_HIGH2LOW_TRANSFORM2 transforms a high .sel into a low .sel (or
%  inverse) by using a .htl file to find the corresponding high or low
%  particles
%
%
%PARAMETERS
%
%  INPUT
%   input_sel            xmipp sel filename
%   input_doc            (opt.) xmipp doc filename
%   input_htl            htl filename
%   output_sel           filename of the output .sel containig the tranfered
%   f_output_doc         (opt.) filenam of output doc-file
%   dir_flag             high2low,low2high
%   all_flag             (0) get high and low instead of the correspinding high or low
%
%
%  OUTPUT
%
%EXAMPLE
%
%
% tom_av2_high2low_transform2('15_corrf_high_128.sel','','15_corrf_high_128_clean.mat_15_corrf_low_128.mat.htl','test_low.sel','','high2low');
% %converts 15_corrf_high_128.sel 2 test_low.sel direction is high2low
%
% tom_av2_high2low_transform2('15_corrf_high_128.sel','in_doc_high.doc','15_corrf_high_128_clean.mat_15_corrf_low_128.mat.htl','test_low.sel','out_doc.doc','high2low');
% %converts 15_corrf_high_128.sel 2 test_low.sel direction is high2low and
% %replaces particle names in the doc file ...to use for rough estimation
% % of tha angles with the high images
%
%
% tom_av2_high2low_transform2('15_corrf_high_128.sel','in_doc_high.doc','15_corrf_high_128_clean.mat_15_corrf_low_128.mat.htl','test_low.sel','out_doc.doc','high2low',1);
% %get high and low instead of the correspinding high or low ...using high and low in one refinement ...bootstrap!
%
% NOTE:
%
% input_htl is a 4 column textfile with path_high path_low part_idx_high part_idx_low
%
% tail 16_corrf_high_128.mat_16_corrf_low_128.mat.htl
%
% /fs/pool/pool-nickell/26S/em/data/bohn/2d/100119_p47f11/log/rec/parts_high_128/parts_609090.spi /fs/pool/pool-nickell/26S/em/data/bohn/2d/100119_p47f11/log/rec/parts_low_128/parts_631378.spi 609090 631378
% /fs/pool/pool-nickell/26S/em/data/bohn/2d/100119_p47f11/log/rec/parts_high_128/parts_609091.spi /fs/pool/pool-nickell/26S/em/data/bohn/2d/100119_p47f11/log/rec/parts_low_128/parts_631379.spi 609091 631379
% /fs/pool/pool-nickell/26S/em/data/bohn/2d/100119_p47f11/log/rec/parts_high_128/parts_609092.spi /fs/pool/pool-nickell/26S/em/data/bohn/2d/100119_p47f11/log/rec/parts_low_128/parts_631380.spi 609092 631380
% /fs/pool/pool-nickell/26S/em/data/bohn/2d/100119_p47f11/log/rec/parts_high_128/parts_609093.spi /fs/pool/pool-nickell/26S/em/data/bohn/2d/100119_p47f11/log/rec/parts_low_128/parts_631381.spi 609093 631381
% /fs/pool/pool-nickell/26S/em/data/bohn/2d/100119_p47f11/log/rec/parts_high_128/parts_609094.spi /fs/pool/pool-nickell/26S/em/data/bohn/2d/100119_p47f11/log/rec/parts_low_128/parts_631382.spi 609094 631382
% /fs/pool/pool-nickell/26S/em/data/bohn/2d/100119_p47f11/log/rec/parts_high_128/parts_609095.spi /fs/pool/pool-nickell/26S/em/data/bohn/2d/100119_p47f11/log/rec/parts_low_128/parts_631383.spi 609095 631383
% /fs/pool/pool-nickell/26S/em/data/bohn/2d/100119_p47f11/log/rec/parts_high_128/parts_609096.spi /fs/pool/pool-nickell/26S/em/data/bohn/2d/100119_p47f11/log/rec/parts_low_128/parts_631384.spi 609096 631384
% /fs/pool/pool-nickell/26S/em/data/bohn/2d/100119_p47f11/log/rec/parts_high_128/parts_609097.spi /fs/pool/pool-nickell/26S/em/data/bohn/2d/100119_p47f11/log/rec/parts_low_128/parts_631385.spi 609097 631385
% /fs/pool/pool-nickell/26S/em/data/bohn/2d/100119_p47f11/log/rec/parts_high_128/parts_609098.spi /fs/pool/pool-nickell/26S/em/data/bohn/2d/100119_p47f11/log/rec/parts_low_128/parts_631386.spi 609098 631386
%
%
%
%
%REFERENCES
%
%SEE ALSO
%
%  tom_av2_high2low_match2.m
%
%   created by  fb
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

if (nargin <6)
    error('dir flag not defined!!');
end;


if (nargin <7)
    all_flag=0;
end;

if (isempty(f_output_doc)==0)
    unix(['head -1 ' input_doc ' > ' f_output_doc]);
end;


disp(' ');
disp(['reading ' input_sel]);
in_sel=importdata(input_sel);
disp(['reading ' input_htl]);
in_htl=importdata(input_htl);

if (isempty(input_doc)==0)
    disp(['reading doc...' input_doc]);
    doc=tom_xmippdocread(input_doc);
end

disp(' ');


i=length(in_htl.data(:,1));
names_h=in_htl.textdata(:,1);
names_l=in_htl.textdata(:,2);
names_h_tmp=in_htl.data(:,1);
names_l_tmp=in_htl.data(:,2);

disp(' ');

if (length(unique(names_l_tmp) ) ==i)
    disp('Low column part numbers in htl are unique !');
else
    disp(['length uniuqe: ' length(unique(names_l_tmp)) ' length totoal: ' length(names_l_tmp)]);
    [isunique lines]=tom_av2_xmipp_check_unique(input_htl);
    error('Low column part numbers in htl are not unique !');
end;

if (length(unique(names_h_tmp) ) ==i)
    disp('High column part numbers in htl are unique !');
else
    [isunique lines]=tom_av2_xmipp_check_unique(input_htl);
    error('high column part numbers in htl are not unique !');
end;

%converting input sel 2 index
in_sel_tmp=zeros(length(in_sel.textdata),1);
for ii=1:length(in_sel.textdata)
    [a b c]=fileparts(in_sel.textdata{ii});
    [a d]=strtok(b,'_');
    in_sel_tmp(ii)=str2double(strrep(d,'_',''));
end;

disp(' ');
if (length(unique(in_sel_tmp) ) ==ii)
    disp('part numbers in .sel are unique !');
else
    %error('part numbers in sel are not unique !');
    disp('part numbers in .sel are NOT unique !');
end;

disp(' ');


in_count=1;
out_count=1;
%build up lookuptable
five_p=round(length(in_sel.textdata)./20);
zz_p=1;
tic;

fp=fopen(output_sel,'w');
all_idx=zeros(length(in_sel.textdata),1);
zz_new_idx=0;

if (isempty(input_doc)==0)
    same_flag=0;
    if  (sum([doc(:).part_idx]-[in_sel_tmp(:)]')==0)
        same_flag=1;
    end;
    if (all_flag==0)
        new_doc=doc;
    else
        new_doc=doc;
        new_doc=cat(1,new_doc,doc);
    end;
end;

%do the matching
for i=1:length(in_sel.textdata)
    
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
    
    if (isempty(f_output_doc)==0 && isempty(idx)==0 )
        if (same_flag==1)
            if (all_flag==1  && strcmp(names_h{idx},names_l{idx})==0)
                zz_new_idx=zz_new_idx+1;
                new_doc(zz_new_idx)=doc(i);
                new_doc(zz_new_idx).name=strrep(doc(i).name,is_name,corr_name2);
            end;
            zz_new_idx=zz_new_idx+1;
            new_doc(zz_new_idx)=doc(i);
            new_doc(zz_new_idx).name=strrep(doc(i).name,is_name,corr_name);
        else
            zz_new_idx=zz_new_idx+1;
            all_idx(i)=find([doc(:).part_idx]==in_sel_tmp(3),1);
        end;
        
    end;
    
    %Ouput
    if (mod(i,five_p)==0)
        toc;
        disp([num2str(zz_p.*5) '% done'] );
        zz_p=zz_p+1;
        tic;
    end;
    
end;
fclose(fp);


if (isempty(f_output_doc)==0)
    disp(['Writing doc ' f_output_doc]);
    tmp=new_doc(1:zz_new_idx);
    tmp(1).header=new_doc(1).header;
    tmp(1).part_idx_unique=new_doc(1).part_idx_unique;
    tom_xmippdocwrite(f_output_doc,tmp);
end;


disp(['Particles in HTL-File: ' num2str(in_count-1) ' Particles not in HTL-File: ' num2str(out_count-1) ]);


