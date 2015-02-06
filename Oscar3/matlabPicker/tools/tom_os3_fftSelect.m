function [dim ft ift] = tom_os3_fftSelect(img)

dim = length(size(img));

if(dim == 2)
    ft = @fft2;
    ift = @ifft2;
else
    ft = @fftn;
    ift = @ifftn;
end;