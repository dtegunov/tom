function [ output_args ] = gtom_ctf_correct_folder( varargin )

o = struct('pixelsize', [1.35e-10, 1.35e-10, 1.35e-10],...
           'cs', [2e-3, 2e-3, 2e-3],...
           'cc', [2.2e-3, 2.2e-3, 2.2e-3],...
           'voltage', [300e3, 300e3, 300e3],...
           'defocus', [-4e-6, 0, 0.1e-6],...
           'astigmatismangle', [0, pi, pi/36],...
           'astigmatism', [0, 1e-6, 0.05e-6],...
           'amplitude', [0.07, 0.07, 0.07],...
           'bfactor', [0, 0, 0],...
           'decaycohill', [0, 0, 0],...
           'decayspread', [0.8, 0.8, 0.8],...
           'dimensions', 512,...
           'maskinner', 8,...
           'maskouter', 128);

paramlist = varargin;
for i=1:2:length(paramlist)
    pname = lower(paramlist{i});
    pvalue = paramlist{i+1};
    
    switch(pname)
        case 'pixelsize'
            o.pixelsize = pvalue.*1e-10;
        case 'cs'
            o.cs = pvalue.*1e-3;
        case 'cc'
            o.cc = pvalue.*1e-3;
        case 'voltage'
            o.voltage = pvalue.*1e3;
        case 'defocus'
            o.defocus = pvalue.*1e-6;
        case 'astigmatismangle'
            o.astigmatismangle = pvalue./180.*pi;
        case 'astigmatism'
            o.astigmatism = pvalue.*1e-6;
        case 'amplitude'
            o.amplitude = pvalue;
        case 'bfactor'
            o.decayksquared = pvalue;
        case 'decaycohill'
            o.decaycohill = pvalue;
        case 'decayspread'
            o.decayspread = pvalue;
        case 'dimensions'
            o.dimensions = pvalue;
        case 'maskinner'
            o.maskinner = pvalue;
        case 'maskouter'
            o.maskouter = pvalue;
    end;
end;

paramsfit = [o.pixelsize; 
             o.cs; 
             o.cc; 
             o.voltage; 
             o.defocus; 
             o.astigmatismangle; 
             o.astigmatism; 
             o.amplitude; 
             o.bfactor;
             o.decaycohill;
             o.decayspread];
paramsfit = paramsfit';
paramskernel = [o.dimensions;
                o.maskinner;
                o.maskouter];
paramskernel = paramskernel';

if ~exist('_corr', 'dir')
    mkdir('_corr');
end;
files = dir('*.mrc');
for f=1:size(files,1)
    filename = files(f).name;
    if exist(['_corr/' filename],'file')
        continue;
    end;
    
    image = tom_mrcread(filename);
    image = image.Value;
    
    [imagecorr, bestfit] = gtom_ctf_fit_correct(image, paramsfit, paramskernel);    
    tom_mrcwrite(imagecorr,'name',['_corr/' filename],'style','classic');
    save(['_corr/' filename '.mat'], 'bestfit');
    
    fprintf([filename ':\n'...
    '  Pixelsize: ' num2str(bestfit(1)/1e-10) ' A\n'...
    '  Cs: ' num2str(bestfit(2)/1e-3) ' mm\n'...
    '  Cc: ' num2str(bestfit(3)/1e-3) ' mm\n'...
    '  Voltage: ' num2str(bestfit(4)/1e3) ' kV\n'...
    '  Defocus: ' num2str(bestfit(5)/1e-6) ' um\n'...
    '  Astigmatism angle: ' num2str(bestfit(6)/pi*180) ' °\n'...
    '  Astigmatism: ' num2str(bestfit(7)/1e-6) ' um\n'...
    '  Amplitude: ' num2str(bestfit(8)) '\n'...
    '  B factor: ' num2str(bestfit(9)) '\n'...
    '  Decay part. coh. ill: ' num2str(bestfit(10)) '\n'...
    '  Decay energy spread: ' num2str(bestfit(11)) '\n']);
    
end;
            
end

