function tom_av2_permute_bg(doc_filename,f_mask3d,name_postfix,mask_bin_rs,demo)
%TOM_AV2_PERMUTE_BG premutes background pixels
%   
%
%  tom_av2_permute_bg((doc_filename,f_mask3d,name_postfix,demo)
%
%  TOM_AV2_PERMUTE_BG premutes background pixels
%  
%
%PARAMETERS
%
%  INPUT
%   doc_filename      *.doc filename 
%   f_mask3d          3d mask filename
%   name_postfix      ('_perm')extension for new particle folders 
%   mask_bin_rs       (0) bin and rescale mask ...speed up if circular use
%                        for a spherical mask use 'circ'
%   demo_mode         (0)
%
%  OUTPUT
%
%EXAMPLE
%     
%  matlabpool open local 8;
%  tom_av2_permute_bg('all_cut.doc','mask3d.em');
%
%
%REFERENCES
%
%SEE ALSO
%   tom_av2_em_classify3d
%
%   created by FB 08/09/09
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

if (nargin < 3)
    name_postfix='_perm';
end;

if (nargin < 4)
    mask_bin_rs=0;
end;

if (nargin < 5)
    demo=0;
end;

mask3d=tom_emread(f_mask3d);
mask3d=mask3d.Value;

sz_mask=size(mask3d);

doc=tom_xmippdocread(doc_filename);

fprintf('%s ', ['Converting Angles: ' ]); 

num_of_entries=size(doc,1);

ref_nr=zeros(num_of_entries,1);
euler_ang_zxz=zeros(num_of_entries,3);
euler_ang_zyz=zeros(num_of_entries,3);
euler_ang_zxz_proj=zeros(num_of_entries,3);

for i=1:num_of_entries
    ref_nr(i)=doc(i).ref;
    [xx,angles] = tom_eulerconvert_xmipp(doc(i).rot, doc(i).tilt, doc(i).psi);
    euler_ang_zxz(i,:)=angles;
    euler_ang_zyz(i,:)=[doc(i).rot doc(i).tilt doc(i).psi];
    [aa tttemp]=tom_eulerconvert_xmipp(doc(i).rot, doc(i).tilt, 0);
    euler_ang_zxz_proj(i,:)=tttemp;
    if (mod(i,2000)==0)    
        fprintf('%s','.');
    end;
end;
fprintf('%s \n','done!');


if demo==1
    figure;
end;

 %disp('Warning Filter added! ');
warning off;

if (strcmp(mask_bin_rs,'circ')==1)
    mask=squeeze(sum(mask3d,3))>0;
    parfor i=1:length(doc)
        
        try
            
            tmp_im=tom_spiderread(doc(i).name);
            tmp_im=tmp_im.Value;
            
            %rotate Back!
            tmp_sh=round([-doc(i).xoff -doc(i).yoff]);
            if (doc(i).flip==0)
                mask_alg=tom_move(tom_rotate(mask, -euler_ang_zyz(i,3)),tmp_sh);
            else
                mask_alg=(tom_move(tom_rotate(mask,-euler_ang_zyz(i,3)),[-tmp_sh(1) tmp_sh(2)]));
                mask_alg=tom_mirror(mask_alg,'x');
            end;
            
            if demo==1
                SUBPLOT(5,1,1); tom_imagesc(mask_alg);
                SUBPLOT(5,1,2); tom_imagesc(tom_filter(tmp_im,8));
                SUBPLOT(5,1,3); tom_imagesc(tom_filter(tmp_im,8).*mask_alg);
                SUBPLOT(5,1,4); tom_imagesc(tmp_im);
            end;
            
            indd=find(mask_alg < 0.1);
            ind_rand=randperm(length(indd));
            tmp_im(indd)=tmp_im(indd(ind_rand));
            
            if (demo==1)
                SUBPLOT(5,1,5); tom_imagesc(tmp_im);
            end;
            
            [a b c]=fileparts(doc(i).name);
            %new_name=[a '_perm/' b c];
            new_name=[a name_postfix '/' b c];
            if (strcmp(new_name,doc(i).name)==1)
                disp(['overwriting particles: ' doc(i).name ' == ' new_name ])
                error('new name == old name');
            end;
            warning off;
            %mkdir([a '_perm/']);
            mkdir([a name_postfix '/']);
            warning on;
            tom_spiderwrite(new_name,tmp_im);
            
            
        catch ME
            disp(['Error: ' doc(i).name])
            disp(ME.message);
        end;
        if (mod(i,10000)==0)
            disp([num2str(i) ' particles processed!']);
        end;
        
    end;
    
    
else
    parfor i=1:length(doc)
        
        if (strcmp(mask_bin_rs,'circ')==0)
            if (mask_bin_rs==0)
                mask=squeeze(sum(tom_rotate(mask3d,euler_ang_zxz_proj(i,:)),3))>0;
            else
                mask=imresize(squeeze(sum(tom_rotate(tom_bin(mask3d,mask_bin_rs),euler_ang_zxz_proj(i,:)),3))>0,[sz_mask(1) sz_mask(2)])>0;
            end;
        end;
        
        try
            
            tmp_im=tom_spiderread(doc(i).name);
            %filter image
            
           
           % tmp_im=tom_bandpass(tmp_im.Value,2,45,3);
            
            tmp_im=tmp_im.Value;
           
            %rotate Back!
            tmp_sh=round([-doc(i).xoff -doc(i).yoff]);
            if (doc(i).flip==0)
                mask_alg=tom_move(tom_rotate(mask, -euler_ang_zyz(i,3)),tmp_sh);
            else
                mask_alg=(tom_move(tom_rotate(mask,-euler_ang_zyz(i,3)),[-tmp_sh(1) tmp_sh(2)]));
                mask_alg=tom_mirror(mask_alg,'x');
            end;
            
            if demo==1
                SUBPLOT(5,1,1); tom_imagesc(mask_alg);
                SUBPLOT(5,1,2); tom_imagesc(tom_filter(tmp_im,8));
                SUBPLOT(5,1,3); tom_imagesc(tom_filter(tmp_im,8).*mask_alg);
                SUBPLOT(5,1,4); tom_imagesc(tmp_im);
            end;
            
            indd=find(mask_alg < 0.1);
            ind_rand=randperm(length(indd));
            tmp_im(indd)=tmp_im(indd(ind_rand));
            
            if (demo==1)
                SUBPLOT(5,1,5); tom_imagesc(tmp_im);
            end;
            
            [a b c]=fileparts(doc(i).name);
            %new_name=[a '_perm/' b c];
            new_name=[a name_postfix '/' b c];
            if (strcmp(new_name,doc(i).name)==1)
                disp(['overwriting particles: ' doc(i).name ' == ' new_name ])
                error('new name == old name');
            end;
            warning off;
            %mkdir([a '_perm/']);
            mkdir([a name_postfix '/']);
            warning on;
            tom_spiderwrite(new_name,tmp_im);
            
            
        catch ME
            disp(['Error: ' doc(i).name])
            disp(ME.message);
        end;
        if (mod(i,10000)==0)
            disp([num2str(i) ' particles processed!']);
        end;
        
    end;
end;


warning on;



















