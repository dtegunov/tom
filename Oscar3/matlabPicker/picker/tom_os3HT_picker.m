%use imgName and tmplName as a cell array.
%the cell array entries must include the full absolute or relative path to the volumes
function results = tom_os3HT_picker(image,template,optionsName,names)
%tom_os3HT_picker 
%   
%   Performs picking of one image by spitting the image into multiple
%   subimages. The correlation is determined for each of the subimages and
%   the results are merged to form the original correlation filter.
%   The classification is performed on one single mashine afterwards.
%   The results written to disk depend on the option structure provided.
%   
%   tom_os3HT_picker(volume,template,optionsName)
%
%PARAMETERS
%
%  INPUT
%   image       - the search image / volume
%   template    - the pattern searched
%   optionsName - if string options will be loaded from defined path
%   names       - cell of image / volume filename and template filename
%  
%  OUTPUT
%   results     - the resulting correlation map
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

%% read options file for template matching

    if(~exist('optionsName'))   
            error('Please provide either an optionsFileName or an options structure. tom_os3HT_picker(...,options,...)');
    end;
    if(ischar(optionsName))
        try
            options = tom_os3_readOptions(optionsName);
        catch
            error(['Error while reading options file, does it exist? Filename: ' optionsName]);
        end;
    else
        options = optionsName;
    end;    
    
    res.flags = options;
%%  if classification enabled, read training stack and align it        
%   the first image of the training stack must be the reference for the
%   alignment
%   
    if(isfield(options,'classification') && options.classification.enabled)
        
        if(ischar(options.classification.trainingStack))
            trainingStack = tom_emread(options.classification.trainingStack);
            trainingStack = trainingStack.Value;
        end;

        if(~isfield(options.classification,'alignMask'))
            options.classification.alignMask = tom_os3_sphereMask(template);
        end;
        
  %      trainingStack = tom_os3_alignStack(trainingStack,template,options.classification.alignMask);
        options.classification.trainingStack = trainingStack; 
        clear trainingStack;
    end;    

%%  skip distributed classification
    options.classification.enabled = false;
    options.result.returnResults = false;
    options.result.saveResults = true;
    
    
    if(~options.classification.enabled)
        
        optionsORIG = options;
        options.classification = rmfield(options.classification,'trainingStack');
    end;
%%  save image to disk -> allow workers a subregion read out
%    tom_emwrite(names{1},image);
%%  create list of parallel jobs
    jobs = tom_os3_createJobList(options,image,template,names{1},names{2});
%     tic;
    if(strcmp(options.job.mode,'parallel'))
        
        parfor(jobIterator = 1:length(jobs))
            matchingResults{jobIterator} = tom_os3_findTemplate(jobs{jobIterator});
        end;
        
%%  debug mode    
    else
        
        matchingResults = {};
        for jobIterator = 1:length(jobs)
            matchingResults{jobIterator} = tom_os3_findTemplate(jobs{jobIterator});
        end;
        
    end;
        
    
    %merge distributed results here and write them to disk
    if(~ options.result.returnResults)
        element = matchingResults{1};
        matchingResults = element.dir;
    end;
    
    result = tom_os3_merge2dImage(matchingResults,options);
%     toc;
    if(~options.classification.enabled)
        options = optionsORIG;
    end;    
%%  enable classification on one machine

    options.classification.enabled = true;
    job = result.job;
    %process classification
    if(isfield(options,'classification') && options.classification.enabled)

        peaks = options.filter.xcf * result.ccc + options.filter.psr * result.psr + options.filter.soc * result.autoc;
        pickList = tom_os3_returnPicklist(peaks,result.angles,result.job,options.filter.numberParticles);
        particleStack  = tom_os3_generateParticleStack(pickList,result.job,'',image);

        %align image stack
        template = job.template;

        particleStack = tom_os3_alignStack(particleStack,template,options.classification.alignMask);

        %classify particleStack
        goodFlags = tom_os3_iterativeClassifier(particleStack,options.classification.trainingStack,options);

        pickList = tom_os3_reducePicklist(pickList,goodFlags);

        %save pickList
        posOfSlash = strfind(job.volumeFile,'/');
        if(~isempty(posOfSlash))
            volumeName = job.volumeFile(posOfSlash(end)+1:end);
        else
            volumeName = job.volumeFile;
        end;
        if(isfield(options,'result') && options.result.saveList)
            
            if(~exist([options.job.resultDirectory '/pickLists'],'dir'))
                mkdir([options.job.resultDirectory '/pickLists']);
            end;

            align2d = tom_os3_pickList2Align2d(pickList);
            save([options.job.resultDirectory '/pickLists/align-' volumeName '.mat'],'align2d');
        end;
        
        results.pickList = tom_os3_pickList2Align2d(pickList);
        
        %save stack
        if(isfield(options,'result') && options.result.saveStack)
            particleStack  = tom_os3_generateParticleStack(pickList,result.job,'',image);
            if(~exist([options.job.resultDirectory '/particleStack'],'dir'))
                mkdir([options.job.resultDirectory '/particleStack']);
            end;

            tom_emwrite([options.job.resultDirectory '/particleStack/stack-' volumeName '.em'],particleStack);
        end;
    end;    
    