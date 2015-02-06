function psd = tom_psd_enhance(psd,center,take_log,filter_w1,filter_w2,decay_width,mask_w1,mask_w2)

if nargin < 8
    mask_w2 = 0.3;
end

if nargin < 7
    mask_w1 = 0.025;
end

if nargin < 6
    decay_width = 0.02;
end

if nargin < 5
    filter_w2 = 0.3;
end

if nargin < 4
    filter_w1 = 0.02;
end

%center the psd
if center == true
    psd = fftshift(psd);
end

%take log
if take_log == true
    psd = log10(psd);
end
    
%median filter
psd = imfilter(psd,fspecial('average',3));


psdsize = size(psd,1);
%bandpass
psd = tom_bandpass(psd,filter_w1.*psdsize,filter_w2.*psdsize,decay_width.*psdsize);

%mask
mask = tom_sphere([psdsize psdsize],mask_w2.*psdsize);
mask2 = tom_sphere([psdsize psdsize],mask_w1.*psdsize);
mask2 = (mask2==0);

psd = psd.*mask.*mask2;

%normalize under mask
idx = find(mask==1);
avgval = mean2(psd(idx));
stdval = std2(psd(idx));

psd(idx) = (psd(idx)-avgval)./stdval;

%mask again
mask = tom_sphere([psdsize psdsize],mask_w2.*psdsize.*0.9);
mask2 = tom_sphere([psdsize psdsize],mask_w1.*psdsize);
mask2 = (mask2==0);

psd = psd.*mask.*mask2;

psd = fftshift(psd);
psd(1,1)=0;
psd = fftshift(psd);