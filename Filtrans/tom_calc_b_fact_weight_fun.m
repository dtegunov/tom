function w_func=tom_calc_b_fact_weight_fun(sz,bfactor,objectpixelsize,cutoff)
% TOM_CALC_B_FACT_WEIGHT_FUN calculates a b-factor weight function 
%  
%     w_func=tom_calc_b_fact_weight_fun(sz,b_fact,pixs,cut_off)
%  
%  PARAMETERS
%  
%    INPUT
%     sz                   size of weigting function 
%     bfactor              b factor
%     objectpixelsize      pixelsize in Ang
%     cutoff               (ny-freq )cut off in Ang
%                  
%
%    OUTPUT
%     w_func               weighting function  2d or 3d
%  
%  EXAMPLE
%     %2d:
%     w_func=tom_calc_b_fact_weight_fun([128 128],-400,4.42,10);
%     %3d:
%     w_func=tom_calc_b_fact_weight_fun([128 128 128],-400,4.42,10);
%  
%     w_vol=tom_ap
%  
%  NOTE:
%     
%  
%  REFERENCES
%  
%  SEE ALSO
%
%   tom_apply_weight_function
%  
%     Copyright (c) 2004-2007
%     TOM toolbox for Electron Tomography
%     Max-Planck-Institute of Biochemistry
%     Dept. Molecular Structural Biology
%     82152 Martinsried, Germany
%     http://www.biochem.mpg.de/tom

if (nargin<4)
    cutoff=2.*objectpixelsize;
end;

x=(2*objectpixelsize.*sz(1))./([1:1:sz(1)]);
q_square=1./((x.^2));
apply_range_pixel=((2.*objectpixelsize)./cutoff).*sz(1)./2;

decay_restore=exp(bfactor.*(q_square')./4);

if (length(sz)==1) %1d case
    decay_restore_rs=imresize(decay_restore,((sz(1)./2)./size(decay_restore,1)));
    w_func=decay_restore_rs;
    return;
end;

if (sz(2) > 1 && length(sz) < 3)  %2d case
    decay_restore_rs=imresize(decay_restore,((sz(1)./2)./size(decay_restore,1)));
    min_val=min(decay_restore_rs(:));
    min_val=min_val(1);
    min_val=0.00001;
    mask=tom_spheremask(ones(sz(1),sz(2)),apply_range_pixel,0);
    decay_restore_2d = tom_polar2cart(tom_norm(decay_restore_rs,1));
    zero_idx=find(decay_restore_2d==0);
    %decay_restore_2d(zero_idx)=min_val;
    decay_restore_2d=(1./(decay_restore_2d+min_val)).*mask;
    w_func=decay_restore_2d;
    return;
end;

if (sz(2)>1 && length(sz) > 2) %3d case
    decay_restore_rs=imresize(decay_restore,((sz(1)./2)./size(decay_restore,1)));
    min_val=min(decay_restore_rs(:));
    min_val=min_val(1);
    for ii=1:2.*sz(1)
        for jj=1:sz(1)
            decay_restore_3d_sphere(:,ii,jj) = decay_restore_rs;
        end;
    end;
    decay_restore_3d = tom_sph2cart(decay_restore_3d_sphere);
    zero_idx=find(decay_restore_3d==0);
    decay_restore_3d(zero_idx)=min_val;
    mask=tom_spheremask(ones(sz(1),sz(1),sz(1)),apply_range_pixel,0);
    decay_restore_3d=(1./(decay_restore_3d)).*mask;
    w_func=decay_restore_3d;
    return;
end;




