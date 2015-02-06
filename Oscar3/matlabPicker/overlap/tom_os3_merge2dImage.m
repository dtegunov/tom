function results = tom_os3_merge2dImage(matchingResults,options)

if(ischar(matchingResults))
    matchingResults = readResultsFromFile(matchingResults);
end;

subVolumeSize = options.parallel.subVolumeSize;
if(subVolumeSize == 1 || isempty(subVolumeSize))
    results = matchingResults{1};
    return;
end;

result = matchingResults{1};
job =result.job;

imageSize = job.volumeSize / 2^job.options.modifications.binning;
templateSize = job.templateSize/ 2^job.options.modifications.binning;

resXCF = zeros(imageSize);
resPSR = resXCF;
resSOC = resXCF;
resAng = resXCF;

subImgCounter = 1;
offset = templateSize(1)/2;
%%  step through each subimage and paste the results back
for x = 0:subVolumeSize-1
    for y = 0:subVolumeSize-1
        %calculate position in image
        X1 = 1 + x*imageSize(1)/subVolumeSize;
        X2 = x*imageSize(1)/subVolumeSize + imageSize(1)/subVolumeSize;
        
        Y1 = 1 + y*imageSize(2)/subVolumeSize;
        Y2 = y*imageSize(2)/subVolumeSize + imageSize(2)/subVolumeSize;
        
        %set result variables
        result = matchingResults{subImgCounter};
        xcf = result.ccc;
        psr = result.psr;
        soc = result.autoc;
        ang = result.angles;
        
        %calculate position in subImage
        subX1 = 1;
        if(X1 > 1)
            subX1 = 1 + templateSize(1)/2+offset;
        end;
        
        subX2 = size(xcf,1);
        if(X2 < imageSize(1))
            subX2 = subX2 - templateSize(1)/2-offset;
        end;

        subY1 = 1;
        if(Y1 > 1)
            subY1 = 1 + templateSize(2)/2+offset;
        end;
        
        subY2 = size(xcf,2);
        if(Y2 < imageSize(2))
            subY2 = subY2 - templateSize(2)/2-offset;
        end;        
        
        %paste the current result
        resXCF(X1:X2,Y1:Y2) = xcf(subX1:subX2,subY1:subY2);
        resPSR(X1:X2,Y1:Y2) = psr(subX1:subX2,subY1:subY2);
        resSOC(X1:X2,Y1:Y2) = soc(subX1:subX2,subY1:subY2);
        resAng(X1:X2,Y1:Y2) = ang(subX1:subX2,subY1:subY2);        
  
        %step to next subImage
        subImgCounter = subImgCounter +1;
        
    end;
end;


results.ccc       = single(resXCF);
results.psr       = single(resPSR);
results.autoc     = single(resSOC);
results.angles    = single(resAng);
results.job       = job;



if(isfield(options,'result') && options.result.saveResults)
    written = tom_os3_secureSave(result,job,options,-1,false);
end;

function matchingResults = readResultsFromFile(path)
    matchingResults = {};

%%  check if volumeDirectory ends with a slash
    if(strcmp(path,'/'))
        slash = '';
    else
        slash = '/';
    end;
    
%%  load each result file to memory  
    results = dir([path slash '*mat']);
    results = {results.name};

    for i=1:length(results)
        load([path slash results{i}]);
        matchingResults{i} = result;
    end;
    
    tmpResults = matchingResults;
    
    %sort results by their id
    for i=1:length(matchingResults)
        result = matchingResults{i};
        id = result.job.subVolumeNumber;
        tmpResults{id} = result;
    end;
    
    matchingResults = tmpResults;
    clear tmpResults;