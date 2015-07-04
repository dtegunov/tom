function [ params ] = gtom_ctfparams_bracket( varargin )

o = struct('pixelsize', [0 0 0],...
           'cs', [0 0 0],...
           'voltage', [0 0 0],...
           'defocus', [0 0 0],...
           'astigmatismangle', [0 0 0],...
           'astigmatism', [0 0 0],...
           'amplitude', [0 0 0],...
           'bfactor', [0 0 0],...
           'scale', [0 0 0],...
           'phaseshift', [0 0 0]);

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
            