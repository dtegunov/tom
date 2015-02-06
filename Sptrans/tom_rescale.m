function out=tom_rescale(in,newsize,mask)

in_fft = fftshift(fftn(in));

if size(newsize,2)<3
    insize = size(in);
    insize(3) = 1;
    newsize(3) = 1;
else
    insize = size(in);
end;

out = complex(zeros(newsize), zeros(newsize));
outoffset = floor(newsize./2) - floor(insize./2) + 1;
inoffset = floor(insize./2) - floor(newsize./2) + 1;

if (nargin<3)
    mask = ones(size(out));
end;

out(max(outoffset(1), 1):min(outoffset(1)+insize(1)-1, end), max(outoffset(2), 1):min(outoffset(2)+insize(2)-1, end), max(outoffset(3), 1):min(outoffset(3)+insize(3)-1, end)) = ...
    in_fft(max(inoffset(1), 1):min(inoffset(1)+newsize(1)-1, end), max(inoffset(2), 1):min(inoffset(2)+newsize(2)-1, end), max(inoffset(3), 1):min(inoffset(3)+newsize(3)-1, end));

out = out.*mask;
out = ifftshift(out);
out = real(ifftn(out));
out = out./(insize(1)./newsize(1).*insize(2)./newsize(2).*insize(3)./newsize(3));
    
end