function weight = av3_computeweight(ref, pixsize, threshold)
%
%   weight = av3_computeweight(ref, pixsize, threshold)
%
%   REF         reference - dense material is assumed to have HIGH gray
%               value
%   pixsize     pixelsize of reference in Angstroem
%   threshold   gray value threshold
%
%   weight      weight in kDa
%

%compute weight per pixel - all in meter and kDa
rhoinkg = 1.3/(0.1^3);
rho = rhoinkg * 1.6605402 *10^24;% rho in kDa
volelement = (pixsize*(10^(-10)))^3;
weipp = volelement * rho;
indx = find(ref > threshold);
weight = weipp*size(indx,1);