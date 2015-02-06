%%
function jobs = tom_os3_createJobList(options,volume,template,imageName,templateName)
%% input 
% options        - options structure, formerly known as flags
    jobs = {};
    jobCounter = 0;

    if(length(options) > 1)

	for i=1:length(options)
	    j = tom_os3_createJobList(options{1});	
	end;
	jobs{length(jobs)+1:length(jobs)+length(j)} = j;
    end;

%%  init variables
    parallelType = options.job.jobType;
  

    nodeCount   = options.parallel.nodeCount;
    
    if(ischar(nodeCount))
        nodeCount = str2double(nodeCount);
    end;
    
    if(nodeCount <= 0 || isnan(nodeCount))
        nodeCount = 1;
    end;
    
%%  check if volumeDirectory ends with a slash
    if(strcmp(options.job.volumeDirectory(end),'/'))
        slash = '';
    else
        slash = '/';
    end;

%%  if a filefilter structure has been specified, use it to create a list of volumes for template matching
%   if non filefilter exists, simply read out the directory
    if(~strcmp(options.job.filefilter,'none'))

        load(options.job.filefilter);
        volumes = {};
        for i = 1:length(particlepicker.filefilter)
            
            if(particlepicker.filefilter{i})
                volumes{length(volumes)+1} = [options.job.volumeDirectory slash particlepicker.filelist{i} ];
            end;    
        end;
    else
%%  create volumes list
        v = dir([options.job.volumeDirectory slash '*em']);
        names = {v.name};

        for i = 1:length(names)
            volumes{i} = [options.job.volumeDirectory slash names{i} ];
        end;
    end;
    
    
    
%%  create volumes list
   
    
    if(strcmp(options.job.templateDirectory(end),'/'))
        slash = '';
    else
        slash = '/';
    end; 
    t = dir([options.job.templateDirectory slash '*em' ]);
    names = {t.name};
    for i = 1:length(names)
        templates{i} = [options.job.templateDirectory slash names{i} ];
    end;    
    
%%  create jobs for 3d picking 
    if(strcmp(parallelType,'3d'))
    
        volumes = volumes{1};
        v = tom_emreadc(volumes);
        v.size  = size(v.Value);
        v.value = 0;
        v.file = volumes;
        
        volumes = v;
   
        
%%      for each template create a jobList        
        for(templateIterator = 1:length(templates))    

            template.file   = templates{templateIterator};
            template.value  = tom_emreadc(template.file);
            template.value  = template.value.Value;
            template.size   = size(template.value);
            template.value  = 0;

            subVolumeSize = options.parallel.subVolumeSize;

            subVolumeList = tom_os3_split(zeros(volumes.size),zeros(template.size),subVolumeSize,'coordinates');

            scanArray    = tom_os3_angleList(options,3);        
%%          split angle List into equal parts
            if(nodeCount > 1)
                arrayLength = floor(length(scanArray) / nodeCount);

                i=1;
                counter = 1;
                while(i <= length(scanArray))
                    if(i+arrayLength <=length(scanArray))
                        angleArray{counter} = scanArray(i:i+arrayLength);
                    else
                        angleArray{counter} = scanArray(i:end);
                    end;

                    i = i+arrayLength+1;
                    counter = counter +1;
                end;
            else
                angleArray = {scanArray};
            end;

%%          create job structures here
            for(i =1:length(subVolumeList))
                for a = 1:length(angleArray)
                    jobCounter          = jobCounter + 1;

                    job.volume          = 0;
                    job.volumeSize      = volumes.size;
                    job.volumeFile      = volumes.file;
                    job.template        = 0;
                    job.templateFile    = template.file;
                    job.templateSize    = template.size;
                    job.coordinates     = subVolumeList{i};
                    job.angleList       = angleArray{a};
                    if(a == 1)
                        job.angleOffset     = 0;
                    else    
                        job.angleOffset     = length(angleArray{a-1});
                    end;
                    job.angleListOrig   = scanArray;
                    job.dimension       = 3;
                    job.options         = options;
                    job.id              = jobCounter;

                    jobs{jobCounter}    = job;
                end;
            end;
        end;    
    
    
%%  create jobs for 2d stack picking       
    elseif(strcmp(parallelType,'2d'))
        
        for volumeIterator = 1:length(volumes)
            
            volume.file  =  volumes{volumeIterator};
            volume.value =  0;
            
            for templateIterator = 1:length(templates)
                
                template.file   = templates{templateIterator};
                template.value  = tom_emreadc(template.file);
                template.value  = template.value.Value;
                template.size   = size(template.value);
                template.value  = 0;
                
                scanArray    = tom_os3_angleList(options,2);        

                if(nodeCount > 1)
                    arrayLength = floor(length(scanArray) / nodeCount);

                    i=1;
                    counter = 1;
                    while(i <= length(scanArray))
                        if(i+arrayLength <=length(scanArray))
                            angleArray{counter} = scanArray(i:i+arrayLength);
                        else
                            angleArray{counter} = scanArray(i:end);
                        end;

                        i = i+arrayLength+1;
                        counter = counter +1;
                    end;
                else
                    angleArray = {scanArray};
                end;
                
                
                for a = 1:length(angleArray)
                    jobCounter          = jobCounter + 1;

                    job.volume          = 0;
                    job.volumeFile      = volumes{volumeIterator};
                    job.template        = 0;
                    job.templateFile    = template.file;
                    job.templateSize    = template.size;
                    
                    job.coordinates     = 0;
                    job.angleList       = angleArray{a};
                    job.angleOffset     = 0;
                    job.angleListOrig   = scanArray;
                    job.dimension       = 2;
                    job.options         = options;
                    job.id              = jobCounter;
                    
                    jobs{jobCounter}    = job;
                end;

            end;
        end;
    elseif(strcmp(parallelType,'2d-online'))
        subVolumeSize = options.parallel.subVolumeSize;
        if(numel(subVolumeSize )<1)
            subVolumeSize = 1;
        end;

        subImages = tom_os3_split2dImage(volume,template,subVolumeSize,'coordinates');
%%      for each subImage create one particular job        
        for volumeIterator = 1:length(subImages)
            
%%          for each template, in this case only one                            
            for templateIterator = 1:length(templates)
              
                scanArray    = tom_os3_angleList(options,2);        
%%              split angular list (better do not use it)
                if(nodeCount > 1)
                    arrayLength = floor(length(scanArray) / nodeCount);

                    i=1;
                    counter = 1;
                    while(i <= length(scanArray))
                        if(i+arrayLength <=length(scanArray))
                            angleArray{counter} = scanArray(i:i+arrayLength);
                        else
                            angleArray{counter} = scanArray(i:end);
                        end;

                        i = i+arrayLength+1;
                        counter = counter +1;
                    end;
                else
                    angleArray = {scanArray};
                end;
                
                
                for a = 1:length(angleArray)
                    jobCounter          = jobCounter + 1;

                    job.volume          = subImages{volumeIterator};
                    job.volumeFile      = imageName;
                    job.volumeSize      = size(volume);
                    job.template        = template;
                    job.templateFile    = templateName;
                    job.templateSize    = size(template);
                    
                    job.coordinates     = 0;
                    job.angleList       = angleArray{a};
                    job.angleOffset     = 0;
                    job.angleListOrig   = scanArray;
                    job.dimension       = 2;
                    job.options         = options;
                    job.id              = jobCounter;
                    job.subVolumeNumber = volumeIterator;
                    jobs{jobCounter}    = job;
                end;

            end;
        end;
    else
        errormsg('No type of parallelization attached. Stopping progress, NOW!');
        return;
    end;
