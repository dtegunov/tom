function im = tom_focalpair(image1,image2)

image1 = tom_mtfdeconv(image1);
image2 = tom_mtfdeconv(image2);

im = ifft2(fft2(image1)+fft2(image2));