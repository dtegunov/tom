function options = tom_os3_readOptions(optionsFile)

    optionsFileH = fopen(optionsFile);

%%  ignore comments
    string = fgetl(optionsFileH);
    while(strcmp(string(1),'#'))
        string = fgetl(optionsFileH);
    end;
    
%%  job description
    options.job.jobType = strtrim(parseString(string));
    string = fgetl(optionsFileH);
    options.job.volumeDirectory = strtrim(parseString(string));
    string = fgetl(optionsFileH);
    options.job.templateDirectory = strtrim(parseString(string));
    string = fgetl(optionsFileH);
    options.job.resultDirectory = strtrim(parseString(string));
    string = fgetl(optionsFileH);
    options.job.wisdomDir = strtrim(parseString(string));
    string = fgetl(optionsFileH);
    options.job.mode = strtrim(parseString(string));
    string = fgetl(optionsFileH);
    options.job.filefilter = strtrim(parseString(string));
    
%%  ignore comments
    string = fgetl(optionsFileH);
    while(strcmp(string(1),'#'))
        string = fgetl(optionsFileH);
    end;

%%  correlation properties  
    options.correlation.type = strtrim(parseString(string));
    string = fgetl(optionsFileH);
    options.correlation.maskFile = strtrim(parseString(string));
    string = fgetl(optionsFileH);
    options.correlation.angles.start = angleValues(strtrim(parseString(string)));
    string = fgetl(optionsFileH);
    options.correlation.angles.end = angleValues(strtrim(parseString(string)));
    string = fgetl(optionsFileH);
    options.correlation.angles.increment = angleValues(strtrim(parseString(string)));

    %set angle options to double values
    
    
%%  ignore comments
    string = fgetl(optionsFileH);
    while(strcmp(string(1),'#'))
        string = fgetl(optionsFileH);
    end;

%%  psf options
    options.psf.file = strtrim(parseString(string));
    
%%  ignore comments
    string = fgetl(optionsFileH);
    while(strcmp(string(1),'#'))
        string = fgetl(optionsFileH);
    end;

%%  parallel settings
    options.parallel.jobManager = strtrim(parseString(string));
    string = fgetl(optionsFileH);
    options.parallel.jobName = strtrim(parseString(string));
    string = fgetl(optionsFileH);
    options.parallel.nodeCount = strtrim(parseString(string));
    string = fgetl(optionsFileH);
    options.parallel.subVolumeSize = angleValues(strtrim(parseString(string)));
    
%%  ignore comments
    string = fgetl(optionsFileH);
    while(strcmp(string(1),'#'))
        string = fgetl(optionsFileH);
    end;

%%  modifications
%   bandpass
    options.modifications.bandpass.low  = str2double(parseString(string));
    string = fgetl(optionsFileH);
    options.modifications.bandpass.high = str2double(parseString(string));
    string = fgetl(optionsFileH);
    options.modifications.bandpass.smoothing = str2double(parseString(string));
%binning    
    string = fgetl(optionsFileH);
    options.modifications.binning       = str2double(parseString(string));
    
%%  ignore comments
    string = fgetl(optionsFileH);
    while(strcmp(string(1),'#'))
        string = fgetl(optionsFileH);
    end;    
    
%%  analysis mode
    options.analysis.ccc = str2double(strtrim(parseString(string)));
    string = fgetl(optionsFileH);
    options.analysis.psr = str2double(strtrim(parseString(string)));
    string = fgetl(optionsFileH);
    options.analysis.autocorr = str2double(strtrim(parseString(string)));
    string = fgetl(optionsFileH);
    options.analysis.pce = str2double(strtrim(parseString(string)));
    string = fgetl(optionsFileH);
    options.analysis.confidence = str2double(strtrim(parseString(string)));
    string = fgetl(optionsFileH);
    
%%  close file
    fclose(optionsFileH);
    
%%
    if(~checkValues)
        errordlg('Error while parsing options file. See matlab console for more info.','Options not valid!');
    end;
%%  set other option values    
    


%%
function res = checkValues(options)    
    res = true;
%%
function res = parseString(string);

    p1 = strfind(string,':')+1;
    p2 = strfind(string,';')-1;

    res = (string(p1:p2));

%%
function res = angleValues(string)
    
    if(isnan(str2double(string)))

        [phi rest] = (strtok(string));
        [psi rest] = (strtok(rest));
        theta = (strtrim(rest));

        res = [str2double(phi) str2double(psi) str2double(theta)];
    else
        res = str2double(string);
    end;
