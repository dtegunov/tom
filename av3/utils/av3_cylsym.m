function cylim = av3_cylsym(im)
%
% AV3_CYLSIM can be used to symmetrize a volume cylindrically
%
%   cylim = av3_cylsym(im)
%
%   The cylindrical axis is assumed to be the z-axis.
%
% PARAMETERS
%  INPUT        
%   IM      3d volume
%   
%  OUTPUT
%   CYLSIM  cylindrically symmetrized vol
%
% last change 03/31/05 FF - updated docu
%
pol = tom_cart2cyl(im);
avpol = sum(pol,2)/size(pol,2);
for ind=1:size(pol,2)
    pol(:,ind,:) = avpol;
end;
cylim = tom_cyl2cart(pol);

