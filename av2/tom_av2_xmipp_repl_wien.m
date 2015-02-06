function tom_av2_xmipp_repl_wien(ctf_fold,f_fit_st,img_size,pix_size,wiener_const,wiener_dim,dyn_flag,first_max_corr,b_fun,pre_wien_flag)
%TOM_AV2_XMIPP_REPL_WIEN replaces wiener filter in a xmipp CtfGroups folder
%
%   tom_av2_xmipp_repl_wien(ctf_fold,fit_st)
%
%  TOM_AV2_XMIPP_REPL_WIEN 
%  
%
%PARAMETERS
%
%  INPUT
%   ctf_fold        ctf folder name
%   f_fit_st        filename of fit struct
%   img_size        size of the wiener filter (for padding img_size*pad-factor)
%   pix_size        pixelsize in Ang for the resulting wiener filter (for padding pix_size/ pad-factor)          
%   wiener_const    (0.1) 1./SNR            
%   wiener_dim      (2d) dimension of the wiener filter 2d or 3d 
%   dyn_flag        ('const') flag for dynamic adaption of wiener const
%                   according 2 part number
%   first_max_corr  (0) use 1 for correcting before the first max
%   b_fun           (no b-factor applied) or b-factor
%   pre_wien_flag   (1) apply b-factor before wiener filter     
%
%EXAMPLE
%  
% tom_av2_xmipp_repl_wien('CtfGroups','low_1.em.mat',[128 128],2.1,0.01,0);
%
%REFERENCES
%
%  www.wadsworth.org/spider_doc/spider/docs/techs/ctf/ctf.html#ref
%
%SEE ALSO
%   
% tom_calc_b_fact_weight_fun
%
% NOTE:
% 
%  %wiener defined as:
%  wiener_f = ctf_theory ./ ((ctf_theory.*ctf_theory) + wiener_const);
%
%
%   created by fb ...ole !!
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


if (nargin < 5)
    wiener_const=0.1;
end;

if (nargin < 8)
    first_max_corr=0;
end;



if (nargin < 10)
    pre_wien_flag=1;
end;

%read defocus file
tt=importdata([ctf_fold '/ctf_groups.defocus']);
tmp=tom_spiderread([ctf_fold '/ctf_group000001.ctf']);

%read number of particles per group

if (strcmp(dyn_flag,'const')==0)
    parts_per_group=importdata([ctf_fold '/ctf_groups.imgno']);
end;


sz_im=size(tmp.Value);




load(f_fit_st);
Fit=st_out.Fit;

%pix_size=Fit.EM.Objectpixelsize;
pix_size=pix_size.*1e-10;
voltage=Fit.EM.Voltage;
%img_size=sz_im;
sigma=Fit.decay_det;
q_0_r=Fit.decay_part_coh_ill_det;
Cs=Fit.EM.Cs;
Cc=Fit.EM.Cc;
deltaE=Fit.decay_energy_spread_det;
amplitude_contrast=Fit.amplitude_contrast_det;


if (nargin < 9 || isempty(b_fun))
   if (strcmp(wiener_dim,'3d')) 
         b_fun=ones([sz_im sz_im(1)]);
   else
        b_fun=ones(sz_im);
   end;
    
else
    ttmp=max(size(b_fun));
    if(ttmp(1)==1)
        %b_fun=tom_calc_b_fact_weight_fun(sz_im,b_fun,pix_size.*1e10,(pix_size.*1e10).*2+((pix_size.*1e10).*2).*0.1);
        b_fun=tom_calc_b_fact_weight_fun(sz_im,b_fun,pix_size.*1e10,(pix_size.*1e10).*2);
    end;
end;


mid=floor(sz_im./2)+1;



%inner_b_cut=20;
%inner_b_idx=find(tom_spheremask(ones(sz_im),inner_b_cut));
%b_fun(inner_b_idx)=b_fun(mid(1)+inner_b_cut,mid(1));
%b_fun(b_fun<2)=2;

b_fun(b_fun>30)=30;

wiener_in=wiener_const;

for i=1:size(tt.data,1)
      
    if (i<10)
        base=[ctf_fold '/ctf_group00000'];
    end;
    
    if (i>=10 && i <100)
        base=[ctf_fold '/ctf_group0000'];
    end;
    if (i>=100 )
        base=[ctf_fold '/ctf_group000'];
    end;
     
     
     Dz=tt.data(i,2).*1e-10;
     [phase amplitude decay E_k E_i_k E_e_k]=tom_ctf2(Dz,0,0,pix_size,voltage,img_size,Cs,sigma,q_0_r,Cc,deltaE);
     ctf_theory=(sqrt(1-amplitude_contrast.^2).*phase + amplitude_contrast.*amplitude);
     %ctf_theory=E_e_k.*E_i_k.*decay.*(sqrt(1-amplitude_contrast.^2).*phase + amplitude_contrast.*amplitude);
     
     if (strcmp(wiener_dim,'3d'))
          ctf_l=ctf_theory(floor(sz_im(1)./2)+1, floor(sz_im(1)./2)+1:end);  
          ctf_theory=make_3d_decay_correction(ctf_l',2*pix_size,64); 
          
     end;
     
     ctf_theory=abs(ctf_theory);
     
     
     %no correction up to 1. max
     if (first_max_corr==0)
        
         if (strcmp(wiener_dim,'3d'))
           
         else
            
         end;
        
        if (strcmp(wiener_dim,'3d'))
            inner_cut=min(tom_crossing(diff(ctf_theory(mid(1)+1:end,mid(1),mid(1)) )))+1;
            inner_idx=find(tom_spheremask(ones([sz_im sz_im(1)]),inner_cut));
        else
            inner_cut=min(tom_crossing(diff(ctf_theory(mid(1)+1:end,mid(1)) )))+1;
            inner_idx=find(tom_spheremask(ones(sz_im),inner_cut));
        end;
        ctf_theory(inner_idx)=1;
        
     end;
     
     if strcmp(dyn_flag,'max_down')
        wiener_const=wiener_in.* sqrt(max(parts_per_group.data(:,2)) ./ parts_per_group.data(i,2) );
     end;
     
     
     if (pre_wien_flag==1)
         ctf_theory = ctf_theory .* b_fun;
     end;
     
     disp(['Group num: ' num2str(i) ' parts: ' num2str(parts_per_group.data(i,2)) ' wiener const: ' num2str(wiener_const) ]);
     
     wiener_f = ctf_theory ./ ((ctf_theory.*ctf_theory) + wiener_const);

     
     if (pre_wien_flag==0)
         wiener_f = wiener_f .* b_fun;
     end;
     
     
     if (strcmp(wiener_dim,'3d'))
         tom_spiderwrite([base num2str(i) '.wien'],wiener_f);
     else
         tom_spiderwrite([base num2str(i) '.wien'],fftshift(wiener_f));
     end;
     
end;


function decay_restore_3d=make_3d_decay_correction(decay_restore,cutoff,sx)

%decay_restore=imresize(decay_restore,((sx./2)./size(decay_restore,1)));
for ii=1:2.*sx
    for jj=1:sx
        decay_restore_3d_sphere(:,ii,jj) = decay_restore;
    end;
end;
decay_restore_3d = tom_sph2cart(decay_restore_3d_sphere);



%legacy
%  ctf_abs=abs(ctf_theory);
%  idx_small=find(ctf_abs<tiny_val);
% ctf_abs_corr=ctf_abs;
%ctf_abs_corr(idx_small)=tiny_val;
%tiny_val=0.35;