function allConfigs = distcompUserConfig()
%distcompUserConfig Return all the user configurations for the Distributed 
%Computing Toolbox as a cell array of function handles.
%
%   Copy this file to a directory that is higher on your MATLAB path than
%   MATLABROOT/toolbox/distcomp/user and edit so that it accurately reflects
%   your scheduler and how you want to run your jobs.
%
%   The most common way of modifying this file is to identify your scheduler
%   and customize the corresponding subfunction.
%    - If you are using a job manager, customize the 'jobmanager' subfunction
%    - If you are using LSF, customize the 'lsf' subfunction
%    - If you are using CCS, customize the 'ccs' subfunction
%    - If you are using the generic scheduler interface, customize the 'generic'
%    subfunction
%    - If you are using the mpiexec scheduler interface, customize the 'mpiexec'
%    subfunction
%
%   This API is subject to change in a future release.  Please be aware that the
%   configurations you enter in this file may need to be moved to a new API at a
%   later date.

%   Copyright 2005-2006 The MathWorks, Inc.
%   $Revision: 1.1.10.5 $  $Date: 2006/12/06 01:36:34 $

% In addition to modifying the configurations shown below, you can add new
% configurations.  If you want to create a new configuration called 'myconfig',
% add the subfunction myconfig to this file, and add the function handle
% @myconfig to the following cell array.

allConfigs = {@cluster, @local};

function conf = jobmanager()
%JOBMANAGER Return a sample configuration for a job manager.
%   You might want to change the name of the job manager, and assign values to
%   LookupURL and PathDependencies.

    %% Parameters to the findResource command
    % See the help for findResource for more information.
    % This is required for this configuration to be complete:
    conf.findResource.Type = 'jobmanager';
    % This is required for this configuration to be complete:
    conf.findResource.Name = 'cluster02';
    % This is required for this configuration to be complete:
    conf.findResource.LookupURL = 'cluster02';
    conf.findResource.State = 'running';

    %% Job properties 
    % Use the doc command to obtain more information about these properties.
    % For example:  doc PathDependencies
    conf.job.PathDependencies = {}; 
    conf.job.FileDependencies = {};
    % The following job properties are specific to the job manager:
    conf.job.RestartWorker = false;
    conf.job.MaximumNumberOfWorkers = inf;
    conf.job.MinimumNumberOfWorkers = 1;
    conf.job.Timeout = inf;

    %% Parallel job properties 
    % Use the doc command to obtain more information about these properties.
    % For example:  doc PathDependencies
    conf.paralleljob.PathDependencies = {}; 
    conf.paralleljob.FileDependencies = {};
    conf.paralleljob.MaximumNumberOfWorkers = inf;
    conf.paralleljob.MinimumNumberOfWorkers = 1;
    % The following paralleljob properties are specific to the job manager:
    conf.paralleljob.RestartWorker = false;
    conf.paralleljob.Timeout = inf;
    
    %% Task properties
    % Use the doc command to obtain more information about these properties.
    % For example:  doc CaptureCommandWindowOutput
    conf.task.CaptureCommandWindowOutput = false;
    % The following task properties are specific to the job manager:
    conf.task.Timeout = inf;
    conf.task.FinishedFcn = '';

function conf = local()

    conf.findResource.Type = 'local';
