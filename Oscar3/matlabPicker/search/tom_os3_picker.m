%use imgName and tmplName as a cell array.
%the cell array entries must include the full absolute or relative path to the volumes
function result = tom_os3_picker(optionsName)
%tom_os3_picker
%   
%   tom_os3_picker(optionsName)
%
%PARAMETERS
%
%  INPUT
%   img         - the search image / volume
%   template    - 
%   options     - optional. If not given 
%   
%  
%  OUTPUT
%   res         - the resulting correlation map
%   options     - this updated version of the options structure contains
%                 the statistics (mean and std) of the search image and its
%                 calculated fourier transform. 
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

%measure time
%% generate return structure where all the values are saved
%
%correlation peak
res.peaks = {};
%
%angles of best match
res.angles = {};
%
%peak to sidelobe ratio
res.psr = {};
%
%current template - without rotation
res.templates = {};
%
%the imaga itself
res.img = {};
%
%correlation of peaks and autocorrelation
res.autocorr = {};
%
%storage for names of image and template
res.path.image = {};
res.path.template = {};


%% read options file for template matching

    if(~exist('optionsName'))   
            error('Please provide either a optionsfilename or a options structure. tom_os3_readOptions(filename)');
    end;
    if(ischar(optionsName))

        options = tom_os3_readOptions(optionsName);
    else
        options = optionsName;
    end;    

    if(length(options) > 1)
	optionsArray = options;
	options = options{1};
    end;
    
    res.flags = options;
%% counter
    resCounter = 1;
 
%%  determine dimension
    dimension = options.job.jobType;
    dimension = str2double(dimension(1));
%%
    if(dimension == 2)

%%      create list of parallel jobs
	if(exist('optionsArray'))
        jobs = tom_os3_createJobList(optionsArray);
	else
	    jobs = tom_os3_createJobList(options);
	end;
        if(strcmp(options.job.mode,'parallel'))
% %%          create a job on a worker pool
%             parallelJob = tom_os3_prepareParallel(options.parallel.jobManager,'jobmanager','TaskParallel',options.parallel.jobName);
%         
% %%      create parallel tasks
%             for jobCounter = 1:length(jobs)
%                 createTask(parallelJob,@tom_os3_findTemplate,1,{jobs{jobCounter}});
%             end;
% %%      execute parallel jobs
% 	    disp('STARTING PARALLEL EXECUTION.');            
%             submit(parallelJob);
% 
%             waitForState(parallelJob);
%         tic;
% 	    disp('PARALLEL EXECUTION FINISHED.');
%             matchingResults = getAllOutputArguments(parallelJob);
% 
%             %find errors
%             [result_p errorsum]=tom_get_paraell_results(parallelJob);
%             if (errorsum > 0);
%                 errordlg('An error occurded during parallel execution. See matlab console for more information.','Parallel Error');
%                 for errorCounter=1:size(result_p,1)
%                     disp(['Hostname: ' result_p{errorCounter,1} '  Error: ' result_p{errorCounter,2}]);
%                 end;
%             end;
%         destroy(parallelJob);
%             parfor(jobIterator = 1:length(jobs))
%                 matchingResults{jobIterator} = tom_os3_findTemplate(jobs{jobIterator});
%             end;
% %%      
%         else
            matchingResults = {};
            for jobIterator = 1:length(jobs)
                matchingResults{jobIterator} = tom_os3_findTemplate(jobs{jobIterator});
            end;
        end;

       
%%        
    elseif(dimension == 3)
        
%for each volume loop possible here       

%%      create list of parallel jobs
        %jobs = tom_os3_createJobList(options);
        j = load('/fs/sun17/lv01/pool/pool-nickell/tmp/T2/hoffman/resultBIGJOB/restJobs.mat');
        jobs = j.j;
        
        if(strcmp(options.job.mode,'parallel'))
%%      create a job on a worker pool
% %%      execute parallel jobs            
            matchingResults = {};
            disp(['Number of jobs: ' num2str(length(jobs))]);
            parfor(jobIterator = 1:length(jobs))
                matchingResults{jobIterator} = tom_os3_findTemplate(jobs{jobIterator});
                disp(['Job no : ' num2str(jobIterator) 'submitted']);
            end;
% %         else            
%%      test environment (on one machine)        
%             matchingResults = {};
%             for jobIterator = 1:length(jobs)
%                 matchingResults{jobIterator} = tom_os3_findTemplate(jobs{jobIterator});
%             end;
            
        end;
        
%%      collect results        
        
        result      = tom_os3_collectResults(options);
        
        %picklist    = tom_os3_returnPicklist(tom_os3_peakValue(result.peaks,result.psr,result.autoc,options),result.angles,result.job.templateSize,100,options.modifications.binning);
        
        %subVolumes  = tom_os3_collectSubVolumes(picklist,result.job.volumeFile);
        
        %for i=1:length(subVolumes)
        %    tmp = subVolumes{i};
        %    figure;tom_dspcub(tmp.volume);
        %end

    else
        errordlg('Input volume has no valid dimension. Aborting...');
    end;
