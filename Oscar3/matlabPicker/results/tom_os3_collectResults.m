%%
function results = tom_os3_collectResults(options)
%tom_os3_collectResults
%   
% 
% 
%   tom_os3_collectResults(options)
%
%PARAMETERS
%
%  INPUT
%   options         - options file
%  
%  OUTPUT
%   result          - 
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

    resultPath  = options.job.resultDirectory;
    if(strcmp(resultPath(end),'/'))
        slash ='';
    else
        slash ='/';
    end;
    resultPath = [resultPath slash];
    files       = dir([resultPath '*.mat']);
    files       = {files.name};
    
    if(length(files) == 0)
        error('Result directory is empty! Aborting');
    end;
        

%%  create result structure for each result type
if(strcmp(options.job.jobType,'3d'))
    
%%  get list of result files in the resultPath


    mkdir([resultPath 'subVolumes']);
    load([resultPath files{1}]);
    %result is now awailable in memory
    volumeSize     = result.job.volumeSize;
    %result will be deleted from memory
    clear result;
    
%%  create target volumes    
    cccVolume      = ones(volumeSize,'single')*(-1000);
    anglesVolume   = cccVolume;
    psrVolume      = cccVolume;
    autocVolume    = cccVolume;

%     tom_emwritec('resultCCC.em',volumeSize,'new');
%     tom_emwritec('resultAngles.em',volumeSize,'new');
%     tom_emwritec('resultPSR.em',volumeSize,'new');
%     tom_emwritec('resultAUT.em',volumeSize,'new');    
%%  for each result file    

    for resultIterator = 1:length(files)
        
        %load result structure
        load([resultPath  files{resultIterator}]);
        disp(files{resultIterator});
        system([ 'mv ' resultPath  files{resultIterator} ' ' resultPath 'subVolumes']);
        cccMap     = result.ccc;
        anglesMap  = result.angles;
        psrMap     = result.psr;
        autocorr   = result.autoc;
        
        job        = result.job;     
        
%%      calculate subVolume which will be updated        
        if(~isequal(volumeSize,size(result.ccc)) || ~isequal(job.coordinates.shiftVector,[0,0,0]))
            coordinates    = job.coordinates.coordinatesWithExtension;
            coordinates(1) = coordinates(1) + job.coordinates.shiftVector(1);
            coordinates(2) = coordinates(1) + options.parallel.subVolumeSize(1) - 1;

            coordinates(3) = coordinates(3) + job.coordinates.shiftVector(2);
            coordinates(4) = coordinates(3) + options.parallel.subVolumeSize(2) - 1;

            coordinates(5) = coordinates(5) + job.coordinates.shiftVector(3);
            coordinates(6) = coordinates(5) + options.parallel.subVolumeSize(3) - 1;
        else
            coordinates(1) = 1;
            coordinates(2) = volumeSize(1);

            coordinates(3) = 1;
            coordinates(4) = volumeSize(2);

            coordinates(5) = 1;
            coordinates(6) = volumeSize(3);
        end;
%%      get subarea of target volume
        %its better to read a than the long expression
        a = cccVolume(coordinates(1):coordinates(2),coordinates(3):coordinates(4),coordinates(5):coordinates(6));
%         a=tom_emreadc('resultCCC.em','subregion',[coordinates(1) coordinates(3) coordinates(5)],[options.parallel.subVolumeSize(1) - 1 options.parallel.subVolumeSize(2)-1 options.parallel.subVolumeSize(3)-1]);
%         a=a.Value;
        b = anglesVolume(coordinates(1):coordinates(2),coordinates(3):coordinates(4),coordinates(5):coordinates(6));
%         b=tom_emreadc('resultAngles.em','subregion',[coordinates(1) coordinates(3) coordinates(5)],[options.parallel.subVolumeSize(1) - 1 options.parallel.subVolumeSize(2)-1 options.parallel.subVolumeSize(3)-1]);
%         b=b.Value;
        c = psrVolume(coordinates(1):coordinates(2),coordinates(3):coordinates(4),coordinates(5):coordinates(6));
%         c=tom_emreadc('resultPSR.em','subregion',[coordinates(1) coordinates(3) coordinates(5)],[options.parallel.subVolumeSize(1) - 1 options.parallel.subVolumeSize(2)-1 options.parallel.subVolumeSize(3)-1]);
%         c=c.Value;
        d = autocVolume(coordinates(1):coordinates(2),coordinates(3):coordinates(4),coordinates(5):coordinates(6));
%         d=tom_emreadc('resultAUT.em','subregion',[coordinates(1) coordinates(3) coordinates(5)],[options.parallel.subVolumeSize(1) - 1 options.parallel.subVolumeSize(2)-1 options.parallel.subVolumeSize(3)-1]);
%         d=d.Value;
%         
%%      determine the better peak of each subvolume
        [a b c d] = tom_os3_bestPeak(a , cccMap, b , anglesMap + job.angleOffset, c , psrMap , d , autocorr,job.dimension,options);
        
%%      write back updated values        
        cccVolume(coordinates(1):coordinates(2),coordinates(3):coordinates(4),coordinates(5):coordinates(6))    = a;
        anglesVolume(coordinates(1):coordinates(2),coordinates(3):coordinates(4),coordinates(5):coordinates(6)) = b;
        psrVolume(coordinates(1):coordinates(2),coordinates(3):coordinates(4),coordinates(5):coordinates(6))    = c;
        autocVolume(coordinates(1):coordinates(2),coordinates(3):coordinates(4),coordinates(5):coordinates(6))  = d;
        
%         tom_emwritec('resultCCC.em',a,'subregion',[coordinates(1) coordinates(3) coordinates(5)],[options.parallel.subVolumeSize(1) - 1 options.parallel.subVolumeSize(2)-1 options.parallel.subVolumeSize(3)-1]);
%         tom_emwritec('resultAngles.em',b,'subregion',[coordinates(1) coordinates(3) coordinates(5)],[options.parallel.subVolumeSize(1) - 1 options.parallel.subVolumeSize(2)-1 options.parallel.subVolumeSize(3)-1]);
%         tom_emwritec('resultPSR.em',c,'subregion',[coordinates(1) coordinates(3) coordinates(5)],[options.parallel.subVolumeSize(1) - 1 options.parallel.subVolumeSize(2)-1 options.parallel.subVolumeSize(3)-1]);
%         tom_emwritec('resultAUT.em',d,'subregion',[coordinates(1) coordinates(3) coordinates(5)],[options.parallel.subVolumeSize(1) - 1 options.parallel.subVolumeSize(2)-1 options.parallel.subVolumeSize(3)-1]);
        clear result;
    end;
%%  set return value(s)   
%     results.peaks   = tom_emreadc('resultCCC.em');
%     results.peaks   = results.peaks.Value;
%     results.angles  = tom_emreadc('resultAngles.em');
%     results.angles   = results.angles.Value;
%     results.psr     = tom_emreadc('resultPSR.em');
%     results.psr   = results.psr.Value;
%     results.autoc   = tom_emreadc('resultAUT.em');
%     results.autoc   = results.autoc.Value;
    results.ccc = cccVolume;
    tom_emwrite([resultPath 'ccc.em'],cccVolume);
    clear cccVolume;
    results.psr =psrVolume;
    tom_emwrite([resultPath 'psr.em'],psrVolume);
    clear psrVolume;
    results.autoc = autocVolume;
    tom_emwrite([resultPath 'soc.em'],autocVolume);
    clear autocVolume;
    results.angles = anglesVolume;
    tom_emwrite([resultPath 'ang.em'],anglesVolume);
    clear anglesVolume;
	results.options = options;
    results.job     = job;
    
    
    save([resultPath '3dresult.mat'],'results');
    clear results;
    r.options = options;
    r.job = job;
    results = r;
else
%%  2d collect    
    
%%  for each result of an sourceimage     
    volumeList = {};
    for fileIterator = 1:length(files)
    %create distinct list of images 
        file    = files{fileIterator};
        posOfDashes = strfind(file,'-');
        
        volumeName = file(posOfDashes(2)+1:posOfDashes(3)-1);
        
        %find volumeName in list
        listIterator = 0;
        found = false;
        while(~found && listIterator+1 <=length(volumeList))
            listIterator    = listIterator +1;
            listElement     = volumeList{listIterator};
            found           = strcmp(listElement.name,volumeName);
        end;

        %if loop has been skipped or filename has not been found, init
        %listElement variable
        if(~exist('listElement') || ~found)
            listElement.name  = '';
            listElement.files = {};
        end;
        listElement.name = volumeName;
        listElement.files{length(listElement.files)+1} = file;
        if(listIterator == 0 || ~found)
            volumeList{listIterator+1} = listElement;
        else
            volumeList{listIterator} = listElement;
        end;
    end;

%%  for each of the distinct volumes, create a distinct list of templates used for correlation
%   and fill the newFileList with result names
    newFileList = {};
    for volumeIterator = 1:length(volumeList)
        
        list = volumeList{volumeIterator};
        templateList = {};
        for listIterator = 1:length(list.files)
            
            file = list.files{listIterator};
            
            posOfDashes = strfind(file,'-');
            templateName = file(posOfDashes(3)+1:end-4);
            
            
            %find volumeName in list
            listIterator = 0;
            found = false;
            while(~found && listIterator+1 <=length(templateList))
                listIterator    = listIterator +1;
                listElement     = templateList{listIterator};
                found           = strcmp(listElement.name,templateName);
            end;

            
            if(~exist('listElement') || ~found)
                listElement.name  = '';
                listElement.files = {};
            end;
            listElement.name = templateName;
            listElement.files{length(listElement.files)+1} = file;
            
            if(listIterator == 0 || ~found)
                templateList{listIterator+1} = listElement;
            else
               	templateList{listIterator} = listElement;
            end;
        end;

%%      collect results for each template correlation result        
        for listIterator= 1:length(templateList)
            listElement = templateList{listIterator};
            
            cccVolume = 0;
            psrVolume = cccVolume;
            autocorrVolume = cccVolume;
            anglesVolume = cccVolume;
                        
            for templateIterator = 1:length(listElement.files)
                
                file        = listElement.files{templateIterator};
                load([resultPath  file]);
                job = result.job;
                
                if(size(cccVolume,1) == 1)
                    cccVolume = ones(size(result.ccc),'single') * (-1000);
                    psrVolume = cccVolume;
                    autocorrVolume = cccVolume;
                    anglesVolume = cccVolume;
                end;
                
                [cccVolume anglesVolume psrVolume autocorrVolume] = tom_os3_bestPeak(cccVolume , result.ccc , ... 
                                 anglesVolume , result.angles, ...
                                 psrVolume , result.psr  , ...
                                 autocorrVolume , result.autoc, ...
                                 job.dimension);                
                clear result;
            end;
            

            
            result.ccc       = cccVolume;
            result.psr       = psrVolume;
            result.autoc     = autocorrVolume;
            result.angles    = anglesVolume;
            result.job       = job;

            save([ resultPath  'result-' list.name '-' listElement.name '.mat'],'result');
            newFileList{length(newFileList)+1} = [ resultPath  'result-' list.name '-' listElement.name '.mat'];
            
            clear result;
        end;   
    end;
%%
    
    results.files       = newFileList;
    results.options     = options;
    
end;



