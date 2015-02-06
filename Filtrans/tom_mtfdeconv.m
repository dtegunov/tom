function image = tom_mtfdeconv(image,method,mtf,cutoff)
%TOM_MTFDECONV  Corrects a volume for ctf and mtf
%
%   image = tom_mtfdeconv(image,method,mtf,cutoff)
%
%
%PARAMETERS
%
%  INPUT
%   image           input image or volume with header
%   method          'flip': apply ctf flipping
%                   'mtfdeconv': no ctf correction
%                   'wiener': use the matlab deconvwnr deconvolution
%   mtf             (optional) mtf function (1D), see   /fs/pool/pool-bmsan-apps/tom_dev/data/mtfs
%   cutoff          (optional) cutoff frequency in pixel
%
%  OUTPUT
%   image           corrected volume with header
%
%EXAMPLE
%       
%   image = tom_mtfdeconv(image,'flip',mtf,512);
%
%REFERENCES
%
%SEE ALSO
%
%   created by AK & SN 16/01/07
%
%   Nickell et al., 'TOM software toolbox: acquisition and analysis for electron tomography',
%   Journal of Structural Biology, 149 (2005), 227-234.
%
%   Copyright (c) 2004-2007
%   TOM toolbox for Electron Tomography
%   Max-Planck-Institute of Biochemistry
%   Dept. Molecular Structural Biology
%   82152 Martinsried, Germany
%   http://www.biochem.mpg.de/tom

if ~isstruct(image)
    error('Image structure with header is required.');
end

if nargin < 4
    cutoff = size(image.Value,1);
end

if nargin < 2
    method = 'flip';
end

%correct image for outliers
%image.Value=tom_xraycorrect(image.Value);

if ~strcmp(method,'mtfdeconv') && ~strcmp(method,'crowther')
    %create theoretical ctf
    if ~isempty(image.Header.FocusIncrement) && image.Header.FocusIncrement ~= 0
        [ctf_out] = tom_create_ctf(image.Header.FocusIncrement./10000, image.Value, image.Header.Objectpixelsize/10, image.Header.Voltage/1000, image.Header.Cs);
    else
        [ctf_out] = tom_create_ctf(image.Header.Defocus./10000, image.Value, image.Header.Objectpixelsize/10, image.Header.Voltage/1000, image.Header.Cs);
    end
end

%create mtf+ctf

%simple flipping
if strcmp(method,'flip')
    ctf_correct=tom_limit(ctf_out./abs(ctf_out+0.00001),-1,1);
elseif strcmp(method,'wiener2')
    figure;ctf_out_1d = tom_ctf(image.Header.FocusIncrement./10000,image.Header.Objectpixelsize/10,  image.Header.Voltage/1000, size(image.Value,1), image.Header.Cs, 0.02,2.2,0.8,0);close(gcf);
    qzero= tom_ctfzero(image.Header.FocusIncrement./10000, image.Header.Objectpixelsize/10,  image.Header.Voltage/1000, size(image.Value,1), image.Header.Cs);
    ctfmax = ctf_out_1d(1:round(qzero(1)));
    ctfmax = find(ctfmax==max(ctfmax));
    ctf_out_1d(1:round(ctfmax)) = 1;
    c = 5;
    ctf_correct_1d = abs(ctf_out_1d)./(ctf_out_1d+c);
    ctf_correct = tom_polar2cart(ctf_correct_1d(1:end-1)');
elseif strcmp(method,'wiener')  
    figure;ctf_out_1d = tom_ctf(image.Header.FocusIncrement./10000,image.Header.Objectpixelsize/10,  image.Header.Voltage/1000, size(image.Value,1), image.Header.Cs, 0.02,2.2,0.8,0);close(gcf);
    qzero= tom_ctfzero(image.Header.FocusIncrement./10000, image.Header.Objectpixelsize/10,  image.Header.Voltage/1000, size(image.Value,1), image.Header.Cs);
    ctfmax = ctf_out_1d(1:round(qzero(1)));
    ctfmax = find(ctfmax==max(ctfmax));
    ctf_out_1d(1:round(ctfmax)) = 1;
    ctf_correct = tom_polar2cart(ctf_out_1d(1:end-1)');
    
    
%deconvolution with full ctf function
elseif strcmp(method,'deconvsin')
    ctf_correct=ctf_out+1.000001;
%deconvolution with mtf only
elseif strcmp(method,'mtfdeconv')
    ctf_correct = ones(size(image.Value,1),size(image.Value,2));
%amplitude restoration according to equation 5 in Fernandez and Crowther, Ultramicropscopy, 2006
elseif strcmp(method,'crowther')
    figure;ctf_out = tom_ctf(image.Header.FocusIncrement./10000,image.Header.Objectpixelsize/10,  image.Header.Voltage/1000, size(image.Value,1), image.Header.Cs, 0.02,2.2,0.8,0);close(gcf);
    f1 = mtf(1);
    f2 = mtf(2);
    f1a = round(mtf(1) .* size(image.Value,1)./2);
    f2a = round(mtf(2) .* size(image.Value,1)./2);
    ctf_correct = ones(1,size(image.Value,1)./2);
    
    %region 1
    ctf_correct(1:f2a) = abs(ctf_out(1:f2a))./(ctf_out(1:f2a).^2+f2);
    %region 2
    w = (abs(ctf_out(f2a+1:f1a)) - f2) ./ (f1 - f2);
    ctf_correct(f2a+1:f1a) = w .* (abs(ctf_out(f2a+1:f1a))./f1^2) + (1-w) .* (abs(ctf_out(f2a+1:f1a)) ./ (ctf_out(f2a+1:f1a).^2 - f2));    
    %region 3
    ctf_correct(f1a+1:end) = 1./abs(ctf_out(f1a+1:end-1)+0.0001);
    ctf_correct = tom_polar2cart(ctf_correct');
    
else
    error('Method not implemented.');
end

%implement mtf if given
if nargin > 2 && ~isempty(mtf)
    mtf_img = tom_polar2cart(tom_norm(mtf,1));
    mtf_img = imresize(mtf_img,[size(ctf_correct,1) size(ctf_correct,2)],'bicubic');
    if ~strcmp(method,'wiener')
        ctf_correct = ctf_correct./(mtf_img+.00001);
    else
        ctf_correct = mtf_img.*ctf_correct;
    end
end

ctf_correct = tom_spheremask(ctf_correct,cutoff,cutoff./32);

%create 3d deconvolution function if input is 3d
if size(image.Value,3) > 1
    sx=size(image.Value,1);
    mtf = zeros((sx./2),2.*sx,sx);
    for ii=1:2.*sx
        for jj=1:sx
            mtf(:,ii,jj) = sum(tom_polar2cart(tom_limit(ctf_correct,-1,1)),2); 
        end
    end;
    ctf_correct = tom_sph2cart(mtf);
end

%deconvolute
if size(image.Value,3) > 1
    image.Value=tom_smooth(image.Value,32,10,'zero');
else
    image.Value=tom_smooth(image.Value,32);
end

if strcmp(method,'wiener')
    %nhood = [3 3];
    % Estimate the local mean of f.
    %localMean = filter2(ones(nhood), image.Value) / prod(nhood);
    % Estimate of the local variance of f.
    %localVar = filter2(ones(nhood), image.Value.^2) / prod(nhood) - localMean.^2;
    % Estimate the noise power if necessary.
    %noise = mean2(localVar);

    %image.Value = deconvwnr(image.Value,real(fftshift(ifft2(fftshift(ctf_correct)))),noise);
    image.Value = deconvwnr(image.Value,real(fftshift(ifft2(fftshift(ctf_correct)))));
    image.Value = tom_bandpass(image.Value,0,cutoff,5);
else
    image.Value = tom_apply_weight_function(image.Value,ctf_correct);
end