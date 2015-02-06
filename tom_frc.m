function [ correlation ] = tom_frc( image1, image2 )

imageft = fftshift(fftn(image1));
% imageft = complex(real(imageft),rand(size(imageft)));
ftpolar1 = complex(tom_cart2polar(real(imageft)), tom_cart2polar(imag(imageft)));
imageft = fftshift(fftn(image2));
% imageft = complex(real(imageft),rand(size(imageft)));
ftpolar2 = complex(tom_cart2polar(real(imageft)), tom_cart2polar(imag(imageft)));

num = sum(real(ftpolar1.*conj(ftpolar2)), 2);
denom = sqrt(sum(abs(ftpolar1).^2, 2).*sum(abs(ftpolar2).^2, 2));

correlation = num./denom;

end

