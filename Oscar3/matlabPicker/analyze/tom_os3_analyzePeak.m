function res = tom_os3_analyzePeak(img,peakMap,template,options)
%tom_os3_analyzePeak
%   
%   tom_os3_analyzePeak 
%   
% 
%   tom_os3_analyzePeak(img,peakMap,template,options)
%
%PARAMETERS
%
%  INPUT
%   img         - the search image / volume
%   peakMap     - peakMap for analysis
%   options     - 
%   
%  
%  OUTPUT
%   res         - the resulting correlation map
%
%EXAMPLE
%   tom_amira_createisosurface(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by TH2 07/07/07
%   updated by ..
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

peakMapSize = size(peakMap);
templateSize = size(template);

res = zeros(peakMapSize);
m = tom_os3_max(peakMap);

innerVolumeSize = 4;
sidelobeSize = floor( size(template,1)/2 );
templateVolume = size(template,1)*size(template,2)*size(template,3);
[dim ft ift] = tom_os3_fftSelect(img);

peakMapErased = tom_limit(peakMap,0.8*m,m,'az');

%%
%%calculate Peak To Sidelobe Ratio for each remaining peak according to
%Correlation Pattern Recognition p.150
%PSR calculation is conducted in frequency domain according to Roseman
%fourier space calculation of the correlation coefficent denominator in
%fourier space
if(strcmp(options.analysis.type,'PSR'))

    means = zeros(peakMapSize,'single');
    stds = means;
    res = means;
%     mask = tom_os3_pasteCenter(means,ones(templateSize,'single'));
    mask = tom_os3_pasteCenter(means,tom_os3_sphereMask(template));
    
    if(isfield(options.analysis,'psrSize'))
        zeroRegion = options.analysis.psrSize;
    else
        if(dim == 2)
            zeroRegion  = floor(templateSize/4);
            zeroRegion(3) = 1; 
        elseif(dim == 3)
            zeroRegion = floor(templateSize/4);
        end;
        
        if(min(templateSize/4) > 5)
            if(dim == 2)
                zeroRegion = [5 5 1];
            else
                zeroRegion = [5 5 5];
            end;
        end;
        
        
    end;
    

%     mask = tom_os3_pasteCenter(mask,zeros((zeroRegion)));
    mask = tom_os3_pasteCenter(mask,0*tom_os3_sphereMask(zeros(templateSize/4)));
    means = tom_os3_mean(peakMap,mask);

    stds = tom_os3_std(peakMap,means,mask);
    
    nominator = peakMap - means;
    denominator = stds;

    tol = 1000*eps( max(abs(denominator(:))));
    nonzero = find(denominator > tol);

    res(nonzero) = nominator(nonzero) ./ denominator(nonzero);

%%    
%calculate Peak to Correlation Energy for each peak
%Correlation Pattern Recognition p.150
elseif(strcmp(options.analysis.type,'PCE'))

    fimg = ft(img);
    fpeak = ft(peakMap);
    
    nominator = abs(sum(sum(sum(fimg .* fpeak))))^2;
    denominator = sum(sum(sum(abs(fimg.*fpeak).*abs(fimg.*fpeak))));
    
    res = nominator / denominator;
%%   
%calulate confidence intervals according to fisher transformation
%Penzek - 2005 Image assesment san diego
%http://de.wikipedia.org/wiki/Korrelationskoeffizient#Fisher-Transformation
elseif(strcmp(options.analysis.type,'Confidence'))    

    
    [peakCoor peakValue peakMapErased] = tom_peak(peakMap,sidelobeSize);
    
    z = fisherTrans(peakValue);

    
    %look up q somewhere
    if(options.alpha == 0.05)
        q = 1.96;
    end;
    %calculate left and right value for confidence interval
    left = z - q/sqrt(templateVolume-3);
    right = z + q/sqrt(templateVolume-3);    
    %return interval borders scaled to r value
    res = [peakValue (exp(2*left)-1)/(exp(2*left)+1) (exp(2*right)-1)/(exp(2*right)+1)];
%%
elseif(strcmp(options.analysis.type,'Gradient'))    
    
    
    if(dim == 2)
        [gx,gy] = gradient(peakMap);
        
        gradients = sqrt(gx.*gx + gy.*gy);
        
        peakMask = ones(7,7,1);
    elseif(dim == 3)
        [gx,gy,gz] = gradient(peakMap);
        
        gradients = sqrt(gx.*gx + gy.*gy + gz .* gz);
        peakMask = ones(7,7,7);
     
    end;
    
    
    
    peakMask = tom_os3_pasteCenter(zeros(size(gradients)),peakMask);
    
    smallArea = gradients .* peakMask;
    
    [maxVal maxPos] = tom_os3_max(smallArea);
    
    res.gradients = gradients;
    
    %%
    %get line to plot in the direction of gradient maximum
    %%
    centerX = floor(size(template,1)/2)+1;
    centerY = floor(size(template,2)/2)+1;
    
   
    r = tom_os3_line(peakMap,template,options.analysis.angle);
    res.xLine = r.xLine;
    res.yLine = r.yLine;
    res.zLine = r.zLine;
%%
%calculate correlation between the peakMap calculatet in a prior correlation 
%and the autocorrelation peak of the image itself
elseif(strcmp(options.analysis.type,'Autocorrelation'))

    %unset old values
    options.correlation.calculationAvailable = false;
    %calculate autocorrelation of template
    options.correlation.isAutoc = true;
    if(~isfield(options.correlation,'templateAutocorr'))
        autocorrelation = tom_os3_corr(template,template);
        options.correlation.templateAutocorr = autocorrelation;
    else
        autocorrelation = options.correlation.templateAutocorr;
    end;
    
    %correlate peakMap and autocorrelation peak
    res = tom_os3_corr(peakMap,autocorrelation,options);

end;

%%
%formula for fisher transformation
function res = fisherTrans(r)

res = 0;
%avoid 0s in transformation
if(r == 1)
    return;
end;

res = 0.5*log((1+r)/(1-r));

%%
function res = hasPeak(img)

res = tom_os3_max(img) > 0;
