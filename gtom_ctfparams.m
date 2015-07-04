function [ params ] = gtom_ctfparams( varargin )

o = struct('pixelsize', 1.0e-10,...
           'cs', 2e-3,...
           'voltage', 300e3,...
           'defocus', -3e-6,...
           'astigmatismangle', 0,...
           'astigmatism', 0,...
           'amplitude', 0.07,...
           'bfactor', 0,...
           'scale', 1,...
           'phaseshift', 0.0);

paramlist = varargin;
for i=1:2:length(paramlist)
    pname = lower(paramlist{i});
    pvalue = paramlist{i+1};
    
    switch(pname)
        case 'pixelsize'
            o.pixelsize = pvalue;
        case 'cs'
            o.cs = pvalue;
        case 'voltage'
            o.voltage = pvalue;
        case 'defocus'
            o.defocus = pvalue;
        case 'astigmatismangle'
            o.astigmatismangle = pvalue;
        case 'astigmatism'
            o.astigmatism = pvalue;
        case 'amplitude'
            o.amplitude = pvalue;
        case 'bfactor'
            o.decayksquared = pvalue;
        case 'scale'
            o.scale = pvalue;
        case 'phaseshift'
            o.phaseshift = pvalue;
    end;
end;

params = [o.pixelsize; 
          o.cs; 
          o.voltage; 
          o.defocus; 
          o.astigmatismangle; 
          o.astigmatism; 
          o.amplitude; 
          o.bfactor;
          o.scale;
          o.phaseshift];
            