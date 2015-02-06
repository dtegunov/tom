function [angularCorrelation angularCorrelationVol]= tom_os3_angularCorrelation(job)
%%  init 

    options         = job.options;
    angularCorrelation =[];
    
    if(~isfield(options,'correlation'))
        options.correlation =[];
    end;
%% load volume,template
    if(job.volume == 0)
        volume = tom_emread(job.volumeFile);
        volume = single(volume.Value);
    else
        volume = single(job.volume);
    end;
    
    volumeSize = size(volume);
    
    %volume = tom_norm(volume,'mean0+1std');
    
    if(job.template == 0)
        template = tom_emread(job.templateFile);
        template = single(template.Value);
    else
        template = single(job.template);
    end; 

    
    
    templateSize = size(template);
    
  %  volume = tom_os3_subVolume(position,volume,template,'center');
    
%     if(~isequal(size(volume),templateSize))
%         error('volume and template size are not equal');
%     end;
%%  
    
    %   generate template mask if none has been specified    
    if((~isfield(options.correlation,'mask') || strcmp(options.correlation.maskFile,'none')))
        [ options.correlation.mask options.correlation.maskSize ] = tom_os3_sphereMask(template);
    end;
%   else load the mask    
    if(~isfield(options.correlation,'mask'))
        %determine template name - last strin in job.templateFile
        pos = strfind(job.templateFile,'/');
        tmpName  = job.templateFile(pos(end)+1:end);
        
        if(strcmp(options.correlation.maskFile,'/'))
            slash = '';
        else
            slash = '/';
        end;
        
        mask = tom_emread([options.correlation.maskFile slash 'mask_' tmpName ]);
        options.correlation.mask = mask.Value;
        options.correlation.maskSize = length(find(mask.Value > 0));
    end;
 
    [ tmp options.correlation.templateMean  options.correlation.templateSTD] = tom_os3_normUnderMask(template,options.correlation.mask);


%% determine the center of volume
center = floor(templateSize/2)+1;

if(length(center) == 2)
    center(3) = 1;
end;


cccMap = zeros(templateSize,'single');

angleList = job.angleListOrig;
optionsOrig = options;
angularCorrelationVol = zeros(volumeSize(1),volumeSize(2),'single');
for angleIterator = 1:length(angleList)
            options = optionsOrig;templateSize(1);
            
            if(iscell(angleList))
                angle= angleList{angleIterator};
            else
                angle = angleList(angleIterator);
            end;
        
%%      rotate template and the template mask by angle
            t   = tom_rotate(template,angle,'linear');

%%      apply point spread function and binning to template        
           % [t options]  = tom_os3_modifyImage(t,options);        
%%      calculate ccc peaks, psr map , autocorrelation        
            [cccMap options]    = tom_os3_corr(volume,t,options);
            
            angularCorrelationVol(:,:,angleIterator) = cccMap;
            
%%      get the ccc value in the template center            
            angularCorrelation(length(angularCorrelation)+1) = tom_os3_max(cccMap);
end;