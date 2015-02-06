function im = av3_rotaverage(im)
%
%   im = av3_rotaverage(im)
%
%   AV3_ROTAVERAGE rotationally averages 3D volume IM 
%
spim = tom_cart2sph(im);
rspim = sum(sum(spim,2),3)/size(spim,2)/size(spim,3);
for indphi =1:size(spim,2)
    for indthe =1:size(spim,3)
        spim(:,indphi,indthe)=rspim;
    end;
end;
im = tom_sph2cart(spim);
