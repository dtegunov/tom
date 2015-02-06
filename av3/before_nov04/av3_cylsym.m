function cylim = av3_cylsym(im)

%
%   cylim = av3_cylsym(im)
%
pol = tom_cart2cyl(im);
avpol = sum(pol,2)/size(pol,2);
for ind=1:size(pol,2)
    pol(:,ind,:) = avpol;
end;
cylim = tom_cyl2cart(pol);

