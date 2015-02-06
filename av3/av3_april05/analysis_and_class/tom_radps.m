function radps = tom_radps(im,nsq)
%
%   radps = tom_radps(im,nsq)
%
%   TOM_RADPS computes radially averaged powerspectrum. The input image IM
%   can be subdivided into NSQ squares. The powerspectra of these subimages
%   are added incoherently. 
%
%   FF 10/26/04
if (nargin < 2)
    nsq =1;
end;
dimx = floor(size(im,1)/nsq);ny = floor(size(im,2)/dimx);
for ix=1:nsq
    for iy=1:ny
        x = (ix-1)*dimx;y=(iy-1)*dimx;
        pspart = tom_ps(im(x+1:x+dimx,y+1:y+dimx));
        if (ix ==1 & iy ==1)
            ps = zeros(size(pspart,1),size(pspart,2));
        end;
        ps = ps + pspart;
    end;
end;
pspol = tom_cart2polar(ps/(nsq*ny));
radps = sum(pspol,2)/(size(pspol,2));