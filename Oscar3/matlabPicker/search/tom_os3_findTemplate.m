function results = tom_os3_findTemplate(job)
%tom_os3_findTemplate
%
%   zur zeit heissen alle aufgefuehrten funktionsnamen anders -> muessen noch mit einem NEW
%   versehen werden, sonst gehts nicht
% 
%   tom_os3_findTemplate(job)
%
%   Applies template matching of a template and an given image/volume.
%   The job structure given as the only input parameter specifies the
%   location of the template and the search image/volume.
%   For further information on the job structure see
%   tom_os3_createJobList.
%   
%   FindTemplate iterates on any angle tupel provided by job, rotates the
%   template(tom_rotate), modifies it (tom_os3_modifyImage), correlates the
%   search image/ volume with the template (tom_os3_corr) and writes the
%   results back to the result directory defined in job.
%   The returned results may contain the resulting volumes as well.
%   
%PARAMETERS
%
%  INPUT
%   job     - the job to be executed by tom_os3_findTemplate
%  
%  OUTPUT
%   results - 
%
%EXAMPLE
%   
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   tom_os3_createJobList, tom_os3_corr , tom_os3_
%
%   created by TH2 07/11/07
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

%%  init function
    dimension       = job.dimension;
    angleList       = job.angleList;
    options         = job.options;
    
 
%%  set system wisdom file for fftw        
    hostname = tom_os3_getHostname;

%%  init search volume
    
    if(strcmp(options.job.jobType,'2d-online') && numel(job.volume) == 4)
            
            volume = tom_emreadc3(job.volumeFile,[job.volume(1)-1 , job.volume(3)-1, 0 , (job.volume(2) - job.volume(1))+1, (job.volume(4) - job.volume(3))+1 ,1]);
            volume = single(volume.Value);
            
    elseif(job.volume == 0)
            volume = tom_emreadc3(job.volumeFile);
            volume = single(volume.Value);
    else
        volume = single(job.volume);
    end;
    
 
    if(dimension == 3)
        originalVolumeSize  = size(volume);
    end;
%%  init template    
    if(job.template == 0)
        template = tom_emreadc3(job.templateFile);
        template = single(template.Value);
    else
        template = job.template;
    end;     
%    figure;tom_dspcub(volume);


%%  apply normalization if any exotic correlation is desired

    if(strcmp(options.correlation.type,'MCF'))
        fvolume = fftn(volume);
        try
            fvolume= fvolume ./ sqrt(abs(fvolume));
        catch
            nonzero = find(abs(fvolume) > 0);
            fvolume(nonzero) = fvolume(nonzero) ./ sqrt(abs(fvolume(nonzero)));
        end;

        volume= ifftn(fvolume);
        
        ftemplate = fftn(template);
        
        try
            ftemplate= ftemplate ./ sqrt(abs(ftemplate));
        catch
            nonzero = find(abs(ftemplate) > 0);
            ftemplate(nonzero) = ftemplate(nonzero) ./ sqrt(abs(ftemplate(nonzero)));
        end;
        template = ifftn(ftemplate);
        options.correlation.type = 'FLCF';
        
    elseif(strcmp(options.correlation.type,'POF'))
        fvolume = fftn(volume);
        try
            fvolume= fvolume ./ sqrt(abs(fvolume));
        catch
            nonzero = find(abs(fvolume) > 0);
            fvolume(nonzero) = fvolume(nonzero) ./ (abs(fvolume(nonzero)));
        end;

        volume= ifftn(fvolume);
        
        ftemplate = fftn(template);
        
        try
            ftemplate= ftemplate ./ sqrt(abs(ftemplate));
        catch
            nonzero = find(abs(ftemplate) > 0);
            ftemplate(nonzero) = ftemplate(nonzero) ./ (abs(ftemplate(nonzero)));
        end;
        template = ifftn(ftemplate);
        
        options.correlation.type = 'FLCF';
    end;
%%  select subregion
    if(dimension == 3 && ~isequal(originalVolumeSize,size(template)))


        coordinates     = job.coordinates.coordinatesWithExtension;
        corrdinatesShift= job.coordinates.shiftVector;
        subVolumeSize   = job.coordinates.subVolumeSize;

        volume = volume(coordinates(1):coordinates(2), ...
            coordinates(3):coordinates(4), ...
            coordinates(5):coordinates(6));

        options.correlation.tollerance = 0.00001;
    else
        options.correlation.tollerance = 0.0000000001;
    end;
    if(options.modifications.binning > 0)        
        volume = tom_bin(volume,options.modifications.binning);
    end;
    
    dimension   = tom_os3_fftSelect(volume);
 
    template = tom_norm(100000*(template+1),'mean0+1std');
    
    if(options.modifications.binning > 0)        
        template= tom_bin(template,options.modifications.binning);
    end;
    
%%    
%   generate template mask if none has been specified    
    if((~isfield(options.correlation,'mask') && strcmp(options.correlation.maskFile,'none')))
        [options.correlation.mask options.correlation.maskSize]= tom_os3_sphereMask(template);
    end;
%   else load the mask    
    if(~isfield(options.correlation,'mask'))
        %determine template name 
        pos = strfind(job.templateFile,'/');
        tmpName  = job.templateFile(pos(end)+1:end);
        
        if(strcmp(options.correlation.maskFile,'/'))
            slash = '';
        else
            slash = '/';
        end;
        
        try
            %chech if the template corresponding masks really exists
            mask = tom_emreadc3([options.correlation.maskFile slash 'mask_' tmpName ]);
        catch 
            mask = tom_os3_sphereMask(template);
        end;
        
        if(options.modifications.binning > 0 && ~isequal(size(mask),size(template)))
            mask = tom_bin(mask.Value,options.modifications.binning);
        end;
        
        if(~isequal(size(mask),size(template)))
            error('The template and its mask differ in size. Please check your settings! It might be due to binning, or a wrong template mask.');
        end;
        
        if(isstruct(mask))
            mask = mask.Value;
        end;
        options.correlation.mask = mask;
        options.correlation.maskSize = length(find(mask > 0));
    end;
    
    templateSize = size(template);
    if(length(templateSize) < 3)
        templateSize(3) = 0;
    end;
     


%%  init return values (allocate memory
    if(ischar(options.parallel.subVolumeSize) || ...
       isequal(options.parallel.subVolumeSize,[ 0 0 0]) || ...
       isequal(options.parallel.subVolumeSize,0) || ...
       numel(options.parallel.subVolumeSize) <2 || ...
       options.parallel.subVolumeSize(2) == 0)
        ccc = ones(size(volume)) * (-1000);
    else
        ccc         = ones(options.parallel.subVolumeSize(1), ...
                           options.parallel.subVolumeSize(2), ...    
                           options.parallel.subVolumeSize(3),'single') * (-1000);
    end;
    
    psr         = ccc;
    autocorr    = ccc;
    angles      = ccc;
    if(dimension == 2)
        angularCorrelation = zeros([size(ccc),length(angleList)],'single');
    end;
    
%%  if volume is empty, nothing can be found -> return zeros    
    if(mean(volume(:)) ~= 0 || std(volume(:)) ~= 0)
        
%%      %normalize volume here, -> less calculations in for loop in
          tmp   = tom_os3_modifyImage(template,options);
        [ tmp options.correlation.templateMean  options.correlation.templateSTD] = tom_os3_normUnderMask(tmp,options.correlation.mask);
        
        
%         options.correlation.templateMean = tom_os3_mean(template,options.correlation.mask);
%         options.correlation.templateSTD  = tom_os3_std(template,options.correlation.templateMean,options.correlation.mask);
%         options.correlation.templateMean = options.correlation.templateMean(templateSize(1)/2+1,templateSize(2)/2+1,templateSize(3)/2+1);
%         options.correlation.templateSTD = options.correlation.templateSTD(templateSize(1)/2+1,templateSize(2)/2+1,templateSize(3)/2+1);
%         
        
        options.correlation.imageMean = tom_os3_mean(volume,options.correlation.mask,options.correlation.maskSize);
        options.correlation.imageSTD  = tom_os3_std(volume,options.correlation.imageMean,options.correlation.mask,options.correlation.maskSize);
        options.correlation.fimage    = fftn(volume);  
        options.correlation.calculationAvailable = true;
        

%%  helper variables
        optionsOrig = options;    
        
%%  log
    logFile = fopen([ options.job.resultDirectory '/logFile-' hostname '.txt'],'a');
%%  for each angle do
        for angleIterator = 1:length(angleList)

            if(iscell(angleList))
                angle = angleList{angleIterator};
            else
                angle = angleList(angleIterator);
            end;
        
%%      rotate template and apply template mask
            t   = tom_rotate(template,angle,'linear') .* options.correlation.mask;

%%      apply point spread function and binning to template        
            [t options]  = tom_os3_modifyImage(t,options);
        
%%      calculate ccc peaks, psr map , autocorrelation        
            [cccMap options]    = tom_os3_corr(volume,t,options);

            cccMap = single(cccMap);
            options.analysis.type   = 'PSR';
            psrMap              = single(tom_os3_analyzePeak(volume,cccMap,t,options));

            options.analysis.type   = 'Autocorrelation';

            autocorrMap         = single(tom_os3_analyzePeak(0,cccMap,t,options));
%%      cut out the subregion if any has been specified              
            if(dimension == 3 && ~isequal(originalVolumeSize,size(template)))
                shift       = corrdinatesShift ;

                
                cccMap      = cccMap(1+shift(1):subVolumeSize(1)+shift(1), ...
                                     1+shift(2):subVolumeSize(2)+shift(2), ...
                                     1+shift(3):subVolumeSize(3)+shift(3));
                
                psrMap      = psrMap(1+shift(1):subVolumeSize(1)+shift(1), ...
                                     1+shift(2):subVolumeSize(2)+shift(2), ...
                                     1+shift(3):subVolumeSize(3)+shift(3));
                                 
                autocorrMap = autocorrMap(1+shift(1):subVolumeSize(1)+shift(1), ...
                                          1+shift(2):subVolumeSize(2)+shift(2), ...
                                          1+shift(3):subVolumeSize(3)+shift(3));             
            end;        
            if(dimension == 2)
                angularCorrelation(:,:,angleIterator) = single(cccMap);
            end;
%%        
            [ ccc angles psr autocorr] = tom_os3_bestPeak(ccc,cccMap, ...
                                                          angles,ones(size(angles),'single')*angleIterator, ...
                                                          psr,psrMap, ...
                                                          autocorr,autocorrMap, ...
                                                          dimension,options);
%%      log
            fprintf(logFile,[ num2str(job.id) ' angle: ' mat2str(angle) ' finished\n'  ]);
        end;
%%  if the image/volume is empty, return 0s        
    else
        ccc     = ccc *0;
        psr     = ccc *0;
        autocorr= ccc *0;
    end;

%%  set small values in angularCorrelation to zero
%     if(dimension == 2)
%         angularCorrelation = angularCorrelation .* (abs(angularCorrelation) >0.001);
%     end;


%%  set the results 
%   this job is not returned
    result.ccc       = single(ccc);
    result.psr       = single(psr);
    result.autoc     = single(autocorr);
    result.angles    = single(angles);
    result.job       = job;


%%  apply classification 
    %
    if(isfield(options,'classification')  && options.classification.enabled)

        
        if(ischar(options.classification.trainingStack))
            trainingStack = tom_emreadc3(options.classification.trainingStack);
            trainingStack = trainingStack.Value;
            
            if(~isfield(options.classification,'alignMask'))
                options.classification.alignMask = tom_os3_sphereMask(reference);
            end;
        
            if(isempty(strfind(options.classification.trainingStack,'aligned')))
                trainingStack = tom_os3_alignStack(trainingStack,template,options.classification.alignMask);
            end;
            
            options.classification.trainingStack = trainingStack;
        
        end;
        
        peaks = options.filter.xcf * result.ccc + options.filter.psr * result.psr + options.filter.soc * result.autoc;
        pickList = tom_os3_returnPicklist(peaks,result.angles,result.job,options.filter.numberParticles);
        particleStack  = tom_os3_generateParticleStack(pickList,result.job,'',volume);

        %align image stack
        template = tom_emreadc3(job.templateFile);
        template = single(template.Value);
        particleStack = tom_os3_alignStack(particleStack,template,options.classification.alignMask);

        %classify particleStack
        goodFlags = tom_os3_iterativeClassifier(particleStack,options.classification.trainingStack,options);

        pickList = tom_os3_reducePicklist(pickList,goodFlags);

        %save pickList
        posOfSlash = strfind(job.volumeFile,'/');
        volumeName = job.volumeFile(posOfSlash(end)+1:end);
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
            particleStack  = tom_os3_generateParticleStack(pickList,result.job,'',volume);
            if(~exist([options.job.resultDirectory '/particleStack'],'dir'))
                mkdir([options.job.resultDirectory '/particleStack']);
            end;

            tom_emwrite([options.job.resultDirectory '/particleStack/stack-' volumeName '.em'],particleStack);
        end;
        
        fprintf(logFile,[ num2str(job.id) ' Classification finished\n'  ]);
        
    end;
 
%     job.volume = template;
%     job.template = template;
    
%     if(dimension == 2)    
%         result.angleCorr = single(tom_os3_correlateAngularVolume(angularCorrelation, tom_os3_angularCorrelation(job)));
%     end;

%%  save results if desired
    if((isfield(options,'result') && options.result.saveResults))
        [ written resultFile resultDir ] = tom_os3_secureSave(result,job,options,logFile,true);
    end;
    fclose(logFile);
%%  print to logfile
    logFile = fopen([ options.job.resultDirectory '/logFile.txt'],'a');
    if(~exist('written') || written)
        fprintf(logFile,[ 'Job no ' num2str(job.id) ' is completed.' '\n']);
    else
        fprintf(logFile,[ 'Problems occured when reading result of job no: ' num2str(job.id) '\n' ]);
    end;
    fclose(logFile);
    
%%  this job is returned

    if(isfield(options,'result') && isfield(options.result,'returnResults') && options.result.returnResults)
        results = result;        
    else
        results = job;
        results.status = 1;
%         results.file = resultFile;
%         results.dir  = resultDir;
    end;
    
%%  THE END    
   