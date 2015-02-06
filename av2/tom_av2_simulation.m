function tom_av2_simulation(outputfolder,flag,part_name,model,proj_classes,noise,ctf,mtf,ctf_corr)
% tom_av2_simulation simulates single particle acquisition
%  
%     tom_av2_simulation(outputfolder,flag,model,part_nr,angles,noise,ctf,ctf_corr)
%  
%  tom_av2_simulation simulates single particle acquisition 
%  
%
%
%  PARAMETERS
%  
%    INPUT
%     outputfolder        folder for outputdata
%     flag                output format ('xmipp','eman2',em)
%     part_name           name for file output ('parts_')
%     model               3D-model or cell file with models (sim conf change)
%     proj_classes        angles matrixin zyz !!!!!!!!  
%                         or struct with  proj_classes.doc_name 
%                                         proj_classes.num_of_part
%     noise               (opt) noise struct (default no noise)  
%     ctf                 (opt) ctf-struct (default full transfer (no ctf))  
%     mtf                 (opt) mtf vector in memory 
%     ctf_corr            (opt) ctf_corr struct (default no correction) 
%                          
%    
%   OUTPUT
%
%   model=tom_emread('/fs/scratch/bohn/fsc_calc/simulation/20s_f64.em');
%   model=model.Value;
%   proj_classes=tom_av2_equal_angular_spacing([0 90],[0 359.9],10,'spider');
%   proj_classes(:,3)=0; %no inplane   
%   proj_classes(:,4)=0; %no x-shift    
%   proj_classes(:,5)=0; %no y-shift
%   proj_classes(:,6)=10; % 10 particles per class 
%   proj_classes(:,7)=0; %rand offset per class use 0 to switch off (a little scatchy parameter due to euler angle distribution use 0.4*increment)   
%                        %rand offset is applied to rot and tilt     
%  
%   tom_av2_simulation('./out','xmipp','parts_',model,proj_classes);
%
%   EXAMPLE
%  
%   %20s getting started 
%   %... no inplane rotation no shifts 15 deg incre 10 particles per class 
%   model=tom_emread('/fs/scratch/bohn/fsc_calc/simulation/20s_f64.em');
%   model=model.Value;
%   proj_classes=tom_av2_equal_angular_spacing([0 90],[0 359.9],20,'spider');
%   proj_classes(:,3)=0; %no inplane   -999 for rand
%   proj_classes(:,4)=0; %no x-shift   -999 for rand
%   proj_classes(:,5)=0; %no y-shift   -999 for rand 
%   proj_classes(:,6)=5; % 5 particles per class 
%   proj_classes(:,7)=0; %rand offset per class use 0 to switch off (a little scatchy parameter due to euler angle distribution use 0.4*increment)   
%                        %rand offset is applied to rot and tilt     
%  
%   tom_av2_simulation('./out','xmipp','part_',model,proj_classes);
%
% 
%   %20s advanced with noise model
%   model=tom_emread('/fs/scratch/bohn/fsc_calc/simulation/20s_f64.em');
%   model=model.Value;
%   proj_classes=tom_av2_equal_angular_spacing([0 90],[0 359.9],10,'spider');
%   proj_classes(:,3)=0; %no inplane   -999 for rand  
%   proj_classes(:,4)=0; %no x-shift   -999 for rand   
%   proj_classes(:,5)=0; %no y-shift   -999 for rand 
%   proj_classes(:,6)=10; % 10 particles per class 
%   proj_classes(:,7)=0; %rand offset per class use 0 to switch off (a little scatchy parameter due to euler angle distribution use 0.4*increment)   
%                        %rand offset is applied to rot and tilt     
%   noise.type='gauss';  
%   noise.snr=0.5;  %or use mesured distribution (v) 
% 
%   tom_av2_simulation('./out','xmipp','part_',model,proj_classes,noise); 
%
%
%   %20s advanced with noise model and ctf model
%   model=tom_emread('/fs/scratch/bohn/fsc_calc/simulation/20s_f64.em');
%   model=model.Value;
%   proj_classes=tom_av2_equal_angular_spacing([0 90],[0 359.9],10,'spider');
%   proj_classes(:,3)=0; %no inplane    -999 for rand 
%   proj_classes(:,4)=0; %no x-shift    -999 for rand 
%   proj_classes(:,5)=0; %no y-shift    -999 for rand  
%   proj_classes(:,6)=10; % 10 particles per class 
%   proj_classes(:,7)=0; %rand offset per class use 0 to switch off (a little scatchy parameter due to euler angle distribution use 0.4*increment)   
%                        %rand offset is applied to rot and tilt     
%   noise.type='gauss';  
%   noise.snr=4;  %or use mesured distribution (v) 
%    
%   ctf.Dz=1.0e-6; %in m or distribution v
%   ctf.Dz_Delta=0e-6; %in m or distribution v
%   ctf.Phi_0=0; %angle in deg or distrigbution v
%   ctf.pix_size=0.44e-9; %in m
%   ctf.voltage=300000; %in V
%   ctf.Cs=2e-3; % in m   
%   ctf.amplitude=0.07; %0..1
%
%   tom_av2_simulation('./out','xmipp','part_',model,proj_classes,noise,ctf); 
%
%
%
%
%
%
%
%
%   NOTE
%   angles are used in zyz convention as most popular packages use zyz (em/tom conv is zxz)
%   angles_and_shifts: 
%   rot  tilt  psi  Xoff Yoff parts_per_class rand_fact
%   
%
%
%  REFERENCES
%  
%  SEE ALSO
%     tom_mesure_snr,tom_fit_ctf_gui
%  
%     created by FB 01/24/06
%  
%     Nickell et al., 'TOM software toolbox: acquisition and analysis for electron tomography',
%     Journal of Structural Biology, 149 (2005), 227-234.
%  
%     Copyright (c) 2004-2007
%     TOM toolbox for Electron Tomography
%     Max-Planck-Institute of Biochemistry
%     Dept. Molecular Structural Biology
%     82152 Martinsried, Germany
%     http://www.biochem.mpg.de/tom



if (nargin < 6)
    noise='';
end;

if (nargin < 7)
    ctf='';
end;

if (nargin < 8)
    mtf='';
end;



if (isstruct(proj_classes)==0)
    %read doc file here!
    
end;

%norm model
model=tom_norm(model,'mean0+1std');

if (isempty(ctf)==0)
    if (size(model,1)<64)
        ctf.img_size=[512 512];
    else
        ctf.img_size=[size(model,1) size(model,2)];
    end
    
    ctf.sigma=.2;
    ctf.q_0_r=0.05; 
    ctf.Cc=2.2e-3; 
    ctf.deltaE=0.8; 
end;

if (isempty(mtf)==0)
    sx=size(model,1);
    idx = find(isnan(mtf));  
    mtf(idx) = 0.0001;
    mtf_fct=imresize(mtf,((sx./2)./size(mtf,1)));
    mtf = zeros((sx./2),2.*sx);
    min_mtf=min(mtf_fct);
    for jj=1:sx*2
        mtf(:,jj) = mtf_fct;
    end;
    mtf_cart = tom_polar2cart(mtf);
    %mtf_cart(find(mtf_cart==0))=min_mtf;
    mtf=mtf_cart;
end;



create_file_st(outputfolder,flag);

list=zeros(sum(proj_classes(:,6)),7);

zz=1;
for i=1:size(proj_classes,1) % loop over all angles
   for ii=1:proj_classes(i,6)
        list(zz,:)=[proj_classes(i,1) proj_classes(i,2) proj_classes(i,3) proj_classes(i,4) proj_classes(i,5) zz i];
        zz=zz+1;    
    end;
end;

parfor i=1:size(list,1)
    create_proj(outputfolder,flag,i,model,part_name,list(i,:),noise,ctf,mtf);
    if (mod(i,10)==0)
        disp(['processing ' num2str(i) ]);
    end;
end;

save([outputfolder '/ang_list.mat'],'list');

debug=1;

if (debug==1)
    
    for ii=1:3
        if (ii==1)
            list_tmp=list;
        end;
        if (ii==2)
            list_tmp=list(1:2:end,:);
        end;
        if (ii==3)
            list_tmp=list(2:2:end,:);
        end;
        
        
        parfor i=1:size(proj_classes,1)
            if (ii==1)
                create_classes(outputfolder,list_tmp(i:i+proj_classes(i,6)-1,:),i,ii);
            else
                idx=find(list_tmp(:,7)==i);
                if (isempty(idx)==0)
                    create_classes(outputfolder,list_tmp(idx,:),i,ii);
                end;
            end;
            
        end;
        back_proj(outputfolder,flag,proj_classes,ii);
    end;
    
    
end;


function create_classes(outputfolder,list,class_nr,ii)

im=tom_spiderread([outputfolder '/parts/part_' num2str(list(1,6)) '.spi']);

tmp_sum=im.Value;


for i=2:size(list,1)
    im=tom_spiderread([outputfolder '/parts/part_' num2str(list(i,6)) '.spi']);
    tmp_sum=tmp_sum+im.Value;
end;

if (isempty(i)==1)
    i=1;
end;

tmp_sum=tom_norm((tmp_sum./i)+100,'phase');

if (ii==1)
    tom_spiderwrite([outputfolder '/classes/classes_all/class_' num2str(class_nr) '.spi'],tmp_sum);
end;
 
if (ii==2)
    tom_spiderwrite([outputfolder '/classes/classes_even/class_' num2str(class_nr) '.spi'],tmp_sum);
end;

if (ii==3)
    tom_spiderwrite([outputfolder '/classes/classes_odd/class_' num2str(class_nr) '.spi'],tmp_sum);
end;



function create_proj(outputfolder,flag,nr,model,part_name,param,noise,ctf,mtf)

if (strcmp(flag,'xmipp'))
    
    sz_mod=size(model);
    mask=tom_spheremask(ones(sz_mod),round(sz_mod(1)./2)-1);
    [a tmp_eu] = tom_eulerconvert_xmipp(param(1),param(2),0);
    rot_tmp=sum(tom_rotate(single(tom_norm(model,'mean0+1std').*mask),tmp_eu,'linear'),3);
    
    sz_tmp=size(rot_tmp);
    
    if (param(3)==-999)
        inp_tmp=round((rand(1).*360));
        rot_tmp=tom_rotate(rot_tmp,inp_tmp); 
    end;
    
    if (param(4)==-999 && param(5)==-999)
        sh_tmp=[round((rand(1).*(sz_tmp(1)./21))) round((rand(1).*(sz_tmp(1)./21)))];
        rot_tmp=tom_shift(rot_tmp,sh_tmp); 
    end;
    
    
    
    
    if (isempty(noise)==0)
        if (length(noise.snr)==1)
            snr=noise.snr;
        end;
        rot_tmp=tom_norm(rot_tmp,'mean0+1std');
        var_mod=std2(rot_tmp).^2;
        var_noise=var_mod./snr;
        if (strcmp(noise.type,'gauss'))
           noise_add=tom_norm(randn(size(rot_tmp)),'mean0+1std');
        end;
        if (strcmp(noise.type,'equal'))
           noise_add=tom_norm(rand(size(rot_tmp)),'mean0+1std'); 
        end;
        noise_add=noise_add.*sqrt(var_noise);
        rot_tmp=tom_norm(rot_tmp+noise_add,'mean0+1std');
    end;
    
    if (isempty(ctf)==0 && isempty(mtf))
        [phase amplitude]=tom_ctf2(ctf.Dz,ctf.Dz_Delta,ctf.Phi_0,ctf.pix_size,ctf.voltage,ctf.img_size, ...
         ctf.Cs,ctf.sigma,ctf.q_0_r,ctf.Cc,ctf.deltaE);
%         tmp_im=zeros(ctf.img_size); 
%         pos=round((size(tmp_im)-size(rot_tmp))./2);
%         tmp_im=tom_paste(tmp_im,rot_tmp,pos);
        %mask=tom_spheremask(ones(size(tmp_im)),200,10);
        tmp_im=rot_tmp;
        tmp_im=tom_apply_weight_function(tmp_im,phase+ (ctf.amplitude.*amplitude) );
        rot_tmp=tmp_im;
    end;
    
    
    if (isempty(mtf)==0 && isempty(ctf))
         rot_tmp=tom_apply_weight_function(rot_tmp,mtf);
    end;
    
    if (isempty(ctf)==0 && isempty(mtf)==0)
        [phase amplitude]=tom_ctf2(ctf.Dz,ctf.Dz_Delta,ctf.Phi_0,ctf.pix_size,ctf.voltage,ctf.img_size, ...
         ctf.Cs,ctf.sigma,ctf.q_0_r,ctf.Cc,ctf.deltaE);
%         tmp_im=zeros(ctf.img_size); 
%         pos=round((size(tmp_im)-size(rot_tmp))./2);
%         tmp_im=tom_paste(tmp_im,rot_tmp,pos);
        %mask=tom_spheremask(ones(size(tmp_im)),200,10);
        tmp_im=rot_tmp;
        w_func=(phase +(amplitude.*ctf.amplitude)).*mtf;
        tmp_im=tom_apply_weight_function(tmp_im,w_func);
        rot_tmp=tmp_im;
    end;
    
    
    
    tom_spiderwrite([outputfolder '/parts/' part_name num2str(nr) '.spi'],rot_tmp);
    
    if (nr==1)
        fp=fopen([outputfolder '/part.sel'],'a');
    else
        fp=fopen([outputfolder '/part.sel'],'a');
    end;
    
    fprintf(fp,[outputfolder '/parts/' part_name num2str(nr) '.spi 1\n']);
    fclose(fp);
    
end;


function back_proj(outputfolder,flag,proj_classes,ii)

if (strcmp(flag,'xmipp'))
    
    
    if (ii==1)
        base=[outputfolder '/classes/classes_all/class_'];
    end;
    
    if (ii==2)
        base=[outputfolder '/classes/classes_even/class_'];
    end;
    
    if (ii==3)
        base=[outputfolder '/classes/classes_odd/class_'];
    end;
    
    
    
    for i=1:10000
        try
            im=tom_spiderread([base num2str(i) '.spi']);
            break;
        catch
        end;
    end;
    
    im=im.Value;
    vol=zeros(size(im,1),size(im,1),size(im,1),'single');
    THICK=size(vol,3);
    
    for i=1:size(proj_classes,1)
        [aa tttemp]=tom_eulerconvert_xmipp(proj_classes(i,1),proj_classes(i,2) , 0);
        proj_ang(:,i)=tttemp;
    end;
    
    zz=1;
    tmpp=zeros(3,size(proj_classes,1));
    for i=1:size(proj_classes,1)
        if (exist([base num2str(i) '.spi'],'file'))
            im=tom_spiderread([base num2str(i) '.spi']);
            im=im.Value;
            if (std2(im)~=0 )
                tmpp(:,zz)=[proj_ang(1,i) proj_ang(2,i) proj_ang(3,i)];
                zz=zz+1;
            end;
        end;
    end;
    proj_ang_clean=tmpp(:,1:(zz-1));
    
    
    for i=1:size(proj_ang,2)
        if (exist([base num2str(i) '.spi'],'file'))
            im=tom_spiderread([base num2str(i) '.spi']);
            im=im.Value;
            if (std2(im)~=0 )
                w_func=tom_calc_weight_functionc([size(vol,1) size(vol,2)],proj_ang_clean',THICK,[proj_ang(1,i) proj_ang(2,i) proj_ang(3,i)]);
                w_proj=tom_apply_weight_function(im,w_func);
                tom_backproj3d_euler(vol,w_proj,proj_ang(1,i),proj_ang(2,i),proj_ang(3,i),[0 0 0]);
            end;
        end;
    end;
    
    
     tom_spiderwrite([outputfolder '/models/model_all.spi'],double(vol));
    
     
      
    if (ii==1)
       tom_spiderwrite([outputfolder '/models/model_all.spi'],double(vol));
    end;
    
    if (ii==2)
        tom_spiderwrite([outputfolder '/models/model_even.spi'],double(vol));
    end;
    
    if (ii==3)
        tom_spiderwrite([outputfolder '/models/model_odd.spi'],double(vol));
    end;
     
     
end;




function create_file_st(outputfolder,flag)

if (strcmp(flag,'xmipp'))
    warning off;
    mkdir([outputfolder '/parts/']);
    mkdir([outputfolder '/models/']);
    mkdir([outputfolder '/classes/classes_even/']);
    mkdir([outputfolder '/classes/classes_odd/']);
    mkdir([outputfolder '/classes/classes_all/']);
    warning on;
end;




