function [corrected decay_restore decay_restore_3d]=tom_apply_bfactor(in,objectpixelsize,bfactor,FSC,apply_range)

% [corrected decay_restore decay_restore_3d]=tom_apply_bfactor(in,objectpixelsize,bfactor,FSC,apply_range)
%
x=(2*objectpixelsize.*size(in,1))./([1:1:size(in,1)]);
q_square=1./((x.^2));
if apply_range(2)==Inf
    apply_range_pixel(1)=0;
else
    apply_range_pixel(1)=((2.*objectpixelsize)./apply_range(2)).*size(in,1)./2;
end;
FSC=imresize(FSC,[size(in,1),1]);
FSC_weight=sqrt((2.*FSC)./(1+FSC));
apply_range_pixel(2)=((2.*objectpixelsize)./apply_range(1)).*size(in,1)./2;
decay_restore=FSC_weight.*exp(bfactor.*(q_square')./4);
decay_restore_3d=make_3d_decay_correction(decay_restore,apply_range_pixel,size(in,1));
corrected=tom_apply_weight_function(in,decay_restore_3d);


function decay_restore_3d=make_3d_decay_correction(decay_restore,cutoff,sx)

decay_restore=imresize(decay_restore,((sx./2)./size(decay_restore,1)));
decay_restore_3d_sphere=zeros(length(decay_restore),2.*sx,1*sx); %beck
for ii=1:2.*sx
    for jj=1:sx
        decay_restore_3d_sphere(:,ii,jj) = decay_restore;
    end;
end;
clear('decay_restore');
decay_restore_3d = tom_sph2cart(decay_restore_3d_sphere);
if cutoff(1)==0
    mask_1=0;
else
    mask_1=tom_spheremask(ones(sx,sx,sx),cutoff(1),0);
end;
mask_2=tom_spheremask(ones(sx,sx,sx),cutoff(2),0);
mask=mask_2-mask_1;
clear('mask_2'); %beck
clear('mask_1'); %beck
decay_restore_3d=(1./(decay_restore_3d+0.00001)).*mask;
