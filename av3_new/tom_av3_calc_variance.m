function [vol_var vol_avg vol_var_even vol_var_odd]=tom_av3_calc_variance(basepath,extension,num_of_models,mask,filter_kernel,outputlevel)
%tom_av3_calc_variance calculates a 3d variance Map 
%
%   [vol_var vol_avg vol_var_even vol_var_odd]=tom_av3_calc_variance(basepath,extension,num_of_models,mask,filter_kernel,outputlevel)
%
%   calculates 3d variance of given folder
%
%PARAMETERS
%
%  INPUT
%   basepath           basepath of volumes
%   extension          extension of the volumes  
%   num_of_models      ('all_you_can') number of models
%                       use all_you can to use all models in folder
%   mask               mask 
%   filter_kernel      filter 
%   outputlevel        debug outputlevel (0,1 or 2)
%
%  
%  OUTPUT
%   vol_var           variance volume
%   vol_avg           average volume
%   vol_var_even      even variance volume           
%   vol_var_odd       odd variance volume 
%
%EXAMPLE
%    
%
%  
% [vol_var vol_avg vol_var_even vol_var_odd]=tom_av3_calc_variance({'rec_files4var/vols/rec_'},{'.spi'},3000,'default',2);
%
%
% REFERENCES
%
%   
%
%SEE ALSO
%   tom_xmippsellread, tom_xmippdocread
%
%   created by FB (feat) Heinz Schenk 27/07/07
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


disp('Creating Volume list ...');
zz=1;
for i=1:length(basepath)
    a=fileparts(basepath{i});
    id=strfind(a,'/');
    a_base=a(1:max(id));
    warning off;
    try
    mkdir([a_base '/diffs']);
    catch
    end;
    warning on;
    dd=dir([basepath{i} '*' extension{i}]);
    for ii=1:length(dd)
        f_list{zz}=[a '/' dd(ii).name];
        p_list{zz}=[a_base '/diffs/' dd(ii).name];
        zz=zz+1;
    end;
    all_dd{i}=dd;
end;


for i=1:length(all_dd)
    dd=all_dd{i};
    if (length(dd)==0)
        error(['no volumes found in: ' basepath{i} '*' extension{i} ] );
    end;
end;

spi_flag=1;
if (tom_isemfile(f_list{1}))
    spi_flag=0;
    im_tmp=tom_emreadc(f_list{1});
end;
if (tom_iseman2h5file(f_list{1}))
    spi_flag=2;
    im_tmp=tom_eman2_read(f_list{1});
end;
    
if (spi_flag==1)
    im_tmp=tom_spiderread(f_list{1});
end;

if nargin < 3 
    num_of_models='all_you_can';
end;

if nargin < 4 
    mask=ones(size(im_tmp.Value,1),size(im_tmp.Value,1),size(im_tmp.Value,1));
end;

if (strcmp(mask,'default'))
     mask=tom_spheremask(ones(size(im_tmp.Value,1),size(im_tmp.Value,1),size(im_tmp.Value,1)),round(size(im_tmp.Value,1)./2)-1);
else
    if (isempty(mask))
        mask=ones(size(im_tmp.Value,1),size(im_tmp.Value,1),size(im_tmp.Value,1));
    end;
end;

if nargin < 5
    filter_kernel=0;
end;

if nargin < 6
    outputlevel=1;
end;

if (strcmp(num_of_models,'all_you_can'))
    num_of_models=length(f_list);
end;

if (length(f_list)<num_of_models)
    disp(['found only ' num2str(length(f_list)) ' models' ]);
    disp(['setting num_of_models 2 ' num2str(length(f_list))]);
end;

disp([num2str(num_of_models) ' of ' num2str(length(f_list)) ' models in folder used']);

vol_avg=zeros(size(im_tmp.Value));
disp('Calculting Avg ...');
error_count=0;
five_p=round(num_of_models./20);

zz=1;
tic;



%get statistics
if (num_of_models > 500)
    all_std=zeros(500,1);
    for i=1:500
        vol=my_save_read(f_list{i},spi_flag);
        all_std(i)=std(std(std(vol.Value)));
    end;
    upper=mean(all_std)+(std(all_std).*3.2);
    lower=mean(all_std)-(std(all_std).*3.2);
else
    upper=1000000;
    lower=-1000000;
end;

used_vol_count=0;


new_f_list=f_list;
new_p_list=f_list;
for i=1:num_of_models
    
    try
        vol=my_save_read(f_list{i},spi_flag);
        if (std(std(std(vol.Value))) > upper || std(std(std(vol.Value))) < lower)
            disp([f_list{i} ' skipped bad statistics']);
            continue;
        end;
        vol.Value=double(vol.Value);
        if (filter_kernel~=0)
            vol.Value=tom_filter(vol.Value,filter_kernel);
        end;
        vol.Value=tom_norm(vol.Value,'mean0+1std',mask).*mask;
        
    catch ME
        disp([f_list{i} ' not readable ...skipping!']);
        disp(ME.message);
        error_count=error_count+1;
    end;
    
    if (i>five_p*zz)
        toc;
        disp([num2str(zz.*5) ' % done']);
        tic;
        zz=zz+1;
    end;
    try
        vol_avg=vol_avg+vol.Value;
    catch ME
        disp([f_list{i} ' corrupt ...skipping!']);
        disp(ME.message);
    end;
    
    used_vol_count=used_vol_count+1;
    new_f_list{used_vol_count}=f_list{i};
    new_p_list{used_vol_count}=p_list{i};
    
end;
   
disp([num2str(num_of_models-used_vol_count) ' volumes sorted out!']);


tmp=new_f_list(1:used_vol_count);
clear('new_f_list'); new_f_list=tmp;

tmp=new_p_list(1:used_vol_count);
clear('new_p_list'); new_p_list=tmp;

vol_avg=vol_avg./used_vol_count;

tic;
if (outputlevel==1); fprintf('%s \n','...calculating variance'); end;

[vol_var vol_var_even vol_var_odd]=calc_var(new_f_list,new_p_list,used_vol_count,vol_avg,spi_flag,mask,filter_kernel,outputlevel);

pause(1);

cc=tom_corr(vol_var_even,vol_var_odd,'norm');
[a b]=tom_peak(cc);

pause(1);

if (length(new_f_list) > 3)
    disp(['Odd/Even Variance correlation: ' num2str(b)]);
end;

try
    warning off; mkdir('var_out'); warning on;
    tom_emwrite(['var_out/var_vol.em'],vol_var);
    tom_emwrite(['var_out/var_even.em'],vol_var_even);
    tom_emwrite(['var_out/var_odd.em'],vol_var_odd);
    tom_emwrite(['var_out/avg.em'],vol_avg);
catch ME
    disp('Cannot write variances!');
    disp(ME.message);
end;

toc;


  
        
        
  function [vol_var vol_var_even vol_var_odd]=calc_var(list,p_list,num_of_models,vol_avg,spi_flag,mask,filter_kernel,outputlevel)
%         
%         on_disk=0;
%         try
%             im_big=zeros(size(vol_avg,1),size(vol_avg,1),size(vol_avg,1),num_of_models);
%             disp('Allocation of 4D diff stack successful!');
%         catch ME
%             disp([ME.message]);
%             on_disk=1;
%         end;
%         
%         try
%            parfor i=1:500
%                im_big(:,:,:,i)=zeros(size(mask));
%            end;
%            disp('Serialization 4D diff stack successful!'); 
%         catch ME
%             disp([ME.message]);
%             on_disk=1;
%         end;
        on_disk=0;
        if (num_of_models>100)
            disp('Stack size > 100 calculating var on disk!')
            on_disk=1;
        end;
        
        vol_avg=vol_avg.*mask;
        vol_var=zeros(size(vol_avg));
        vol_var_even=zeros(size(vol_avg));
        vol_var_odd=zeros(size(vol_avg));
        
        idx=1:num_of_models;
        idx_even=find(mod(idx,2).*idx);
        idx_odd=find((mod(idx,2)==0).*idx);
        
        if (outputlevel >=1)
            fprintf('%s ', ['Calculating Variance: ' ]);
        end;
       
        
        if (on_disk==0)
            read_count=zeros(num_of_models,1);
            parfor i=1:num_of_models
                tmp=my_save_read(list{i},spi_flag);
                
                try
                    tmp=double(tmp.Value);
                    if (filter_kernel~=0)
                        tmp=tom_filter(tmp,filter_kernel);
                    end;
                    tmp=tom_norm(tmp,'mean0+1std',mask).*mask;
                    im_big(:,:,:,i)=(tmp-vol_avg).^2;
                    read_count(i)=1;
                catch ME
                    disp([list{i} ' corrupted ...skipping!']);
                    disp([ME.message]);
                end;
                if (outputlevel >=1)
                    fprintf('%s','.');
                end;
            end;
            idx=find(read_count==1);
            idx_even=find(mod(idx,2).*idx);
            idx_odd=find((mod(idx,2)==0).*idx);
            for i=1:length(idx)
                vol_var=vol_var+im_big(:,:,:,idx(i));
            end;
            for i=1:length(idx_even)
                vol_var_even=vol_var_even+im_big(:,:,:,idx_even(i));
            end;
            for i=1:length(idx_odd)
                vol_var_odd=vol_var_odd+im_big(:,:,:,idx_odd(i));
            end;
        
        end;
            
        if (on_disk==1)
        
            five_p=round(num_of_models./20);
            %tic;
            disp(' ');
            disp(['Writing diff Volumes : ' p_list{1}]);
            parfor i=1:num_of_models
                tmp=my_save_read(list{i},spi_flag);
                if (isempty(tmp))
                    continue;
                end;
                tmp=double(tmp.Value);
                if (filter_kernel~=0)
                    tmp=tom_filter(tmp,filter_kernel);
                end;
                tmp=tom_norm(tmp,'mean0+1std',mask).*mask;
                im_var_tmp=(tmp-vol_avg).^2;
                my_save_write(p_list{i},im_var_tmp,0);
                if (mod(i,five_p)==0)
                   % toc;
                    disp([num2str(i) '  done']);
                   % tic;
                end;
            end;
            
            five_p=round(num_of_models./20);
            zz=1;
            diff_error_count=0;
            tic;
            disp(' ');
            disp(['Avg diff Volumes : ' ]);
            for i=1:num_of_models
                im_var_tmp=my_save_read(p_list{i},0);
                if (isempty(im_var_tmp))
                    disp(['skipping ' p_list{i}]);
                    diff_error_count=diff_error_count+1;
                    continue;
                end;
                im_var_tmp=im_var_tmp.Value;
                vol_var=vol_var+im_var_tmp;
                if (find(idx_even==i))
                    vol_var_even=vol_var_even+im_var_tmp;
                end;
                if (find(idx_odd==i))
                    vol_var_odd=vol_var_odd+im_var_tmp;
                end;
                if (i>five_p*zz)
                    toc;
                    disp([num2str(zz.*5) ' % done']);
                    tic;
                    zz=zz+1;
                end;
            end;
            
            disp([num2str(diff_error_count) ' errors in ' num2str(num_of_models) ' diff volumes !']);
        end;
        
        
        
        %norm it norman!
        vol_var=(vol_var./(length(idx)-1));
        vol_var_even=(vol_var_even./(length(idx_even)-1));
        vol_var_odd=(vol_var_odd./(length(idx_odd)-1));
        
        
        
        fprintf('%s \n', ['done! ' ]);
        
        
        function tmp=my_save_read(name,spi_flag)
            
            for i_try=1:15
                try
                    if (spi_flag==0)
                        tmp=tom_emreadc(name);
                    end;
                    if (spi_flag==1)    
                        tmp=tom_spiderread(name);
                    end;
                    if (spi_flag==2)    
                        tmp=tom_eman2_read(name);
                    end;
                    break;
                catch ME
                    disp(['error Reading: ' name ' ' ME.message]);
                    pause2(i_try.*2);
                    tmp='';
                end;
            end;
        
        
            
       function my_save_write(name,vol,spi_flag)
            
           
            for i_try=1:15
                try
                    if (spi_flag==0)
                       tom_emwritec(name,vol);
                    end;
                    if (spi_flag==1)  
                        tom_spiderwrite(name,vol);
                    end;
                    if (spi_flag==2)  
                        tom_eman2_write(name,vol,'volume');
                    end;
                    break;
                catch ME
                    disp(['error writing: ' name ' ' ME.message]);
                    pause2(i_try.*2);
                end;
            end;
        
           
           
           
           
           
           
           
           
           
           
           
           
           
        
        
        