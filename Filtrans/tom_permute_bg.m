function out=tom_permute_bg(in,mask,outputname,grow_rate,num_of_steps,filt_cer,max_error)
%TOM_PERMUTE_BG premutes background pixels
%   
%
%  out=tom_permute_bg(in,mask,outputname,grow_rate,num_of_steps,filt_cer,max_error)
%
%  TOM_PERMUTE_BG premutes background pixels which are outside of a given
%  mask
%  
%
%PARAMETERS
%
%  INPUT
%   in                1d,2d or 3d signal
%   mask              mask to det structure
%   outputname        ('') filename for fileoutput use '' for no hd output
%   grow_rate         (0) grow rate for the initial mask (for smoothing use a rate bigger 1 )
%   num_of_steps      (5) number of growing steps 
%   filt_cer          (3% of mask size) filter for growing the initail mask
%   max_error         (3%) max error for mask growing
%
%  OUTPUT
% 
%  out                permed signal 
%
%EXAMPLE
%     
% im=tom_emread('mona.em');
% mask=tom_spheremask(ones(size(im.Value)),round(min(size(im.Value))./4));
% out=tom_permute_bg(im,mask);
% figure; tom_imagesc(out);
%
% out=tom_permute_bg(im,mask,'tmp.spi',1.008,40,2);
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

in_flag='.em';

if (nargin < 3)
    outputname='';
end;

if (nargin < 4)
    grow_rate=0;
end;

if (nargin < 5)
    num_of_steps=10;
end;

if (nargin < 6)
    filt_cer=10;
end;

if (nargin < 7)
    max_error=10;
end;



if (ischar(in))
    if (tom_isemfile(in))
        in=tom_emread(in);
        in_flag='em';
    else
        in=tom_spiderread(in);
        in_flag='spi';
    end;
end;

if (ischar(mask))
    if (tom_isemfile(mask))
        mask=tom_emread(mask);
    else
        mask=tom_spiderread(mask);
    end;
end;

if (isstruct(in))
    in=in.Value;
end;

if (isstruct(mask))
    mask=mask.Value;
end;

if (grow_rate==0)
    indd=find(mask < 0.1);
    ind_rand=randperm(length(indd));
    in(indd)=in(indd(ind_rand));
else
    
    mask_tmp=mask;
    smooth_ch=(100./num_of_steps):(100./num_of_steps):100;
    std_ch=(4./num_of_steps):(4./num_of_steps):4;
    filt_ch=round(2:-(2./num_of_steps):0);
    for i=1:num_of_steps
        mask_old=mask_tmp;
        mask_tmp=tom_grow_mask(mask_tmp,grow_rate,max_error,filt_cer);
        mask_diff=mask_tmp-mask_old;
        indd=find(mask_diff > 0.1);
        ind_rand=randperm(length(indd));
        ind_rand2=randperm(length(indd));
        cut_len=round(length(ind_rand).*(smooth_ch(i)./100));
        tmp_vox=in(indd(ind_rand(1:cut_len)));
        tmp_vox=clean_stat(tmp_vox,std_ch(i));
        in(indd(ind_rand2(1:cut_len)))=tmp_vox;
        if (filt_ch(i)==0)
            in_f=in;
        else
            in_f=tom_filter(in,filt_ch(i));
        end;
        tmp_vox=in_f(indd(ind_rand(1:cut_len)));
        tmp_vox=clean_stat(tmp_vox,std_ch(i));
        in(indd(ind_rand2(1:cut_len)))=tmp_vox;
        
    end;
    indd=find(mask_tmp < 0.1);
    ind_rand=randperm(length(indd));
    in(indd)=in(indd(ind_rand));
end;


out=in;

if isempty(outputname)==0
    if (strcmp(in_flag,'spi'))
        tom_spiderwrite(outputname,out);
    else
        tom_emwrite(outputname,out);
    end;
end;

        
function out=clean_stat(vox,nstd)

mea_out=mean(vox(:));
std_out=std(vox(:));

idx=find((vox > (mea_out + (nstd.*std_out))) + (vox < (mea_out - (nstd.*std_out))));

out=vox;

out(idx)=mea_out;



        
        
        
        