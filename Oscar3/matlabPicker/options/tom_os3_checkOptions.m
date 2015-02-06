function [valid optionsValid invalidCounter] = tom_os3_checkOptions(options)



invalidCounter =0;
valid = true;

if(strcmp(options.job.jobType,''))
    valid = false;
    options.job.jobType = false;
    invalidCounter =invalidCounter+1;
end;

if(strcmp(options.job.volumeDirectory,''))
    valid = false;
    options.job.jobType = false;
    invalidCounter =invalidCounter+1;
elseif(~exist(options.job.volumeDirectory,'dir'))
    valid = false;
    options.job.volumeDirectory = false;
    invalidCounter =invalidCounter+1;    
end;

if(strcmp(options.job.templateDirectory,''))
    valid = false;
    options.job.templateDirectory = false;
    invalidCounter =invalidCounter+1;
elseif(~exist(options.job.templateDirectory,'dir'))
    valid = false;
    options.job.templateDirectory = false;
    invalidCounter =invalidCounter+1;
end;

if(strcmp(options.job.resultDirectory,''))
    valid = false;
    options.job.resultDirectory = false;
    invalidCounter =invalidCounter+1;
elseif(~exist(options.job.resultDirectory,'dir'))
    valid = false;
    options.job.resultDirectory = false;
    invalidCounter =invalidCounter+1;
end;

if(strcmp(options.job.wisdomDir,''))
    valid = false;
    options.job.wisdomDir = false;
    invalidCounter =invalidCounter+1;
elseif(~exist(options.job.wisdomDir,'dir'))
    valid = false;
    options.job.wisdomDir = false;
    invalidCounter =invalidCounter+1;
end;

if(~strcmp(options.job.filefilter,'none') && ~exist(options.job.filefilter,'file'))
    valid = false;
    options.job.filefilter = false;
    invalidCounter =invalidCounter+1;
end;


optionsValid = options;

% options.correlation.type = 'FLCF';
% options.correlation.maskFile = 'none';
% options.correlation.angles.start =[];
% options.correlation.angles.end = [];
% options.correlation.angles.increment =[];
% 
% options.psf.file = '';
% 
% options.parallel.jobManager ='';
% options.parallel.jobName ='';
% options.parallel.nodeCount ='';
% options.parallel.subVolumeSize =[];
% 
% 
% options.modifications.binning =[];
% options.modifications.bandpass.low = -1;
% options.modifications.bandpass.high = [];
% options.modifications.bandpass.smoothing = [];
% 
% options.analysis.ccc = 1;
% options.analysis.psr = 1;
% options.analysis.autocorr = 1;
% options.analysis.pce = 0;
% options.analysis.confidence = 0;