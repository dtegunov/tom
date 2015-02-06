function cc = tom_ncorr(v1, v2, do_shift)

v1 = (v1 - mean(v1(:))) / std(v1(:), 1);
v2 = (v2 - mean(v2(:))) / std(v2(:), 1);

cc = (real(ifftn(fftn(v1).*conj(fftn(v2))))) / numel(v1);

if (exist('do_shift','var') && do_shift)
    cc = fftshift(cc);
end;


