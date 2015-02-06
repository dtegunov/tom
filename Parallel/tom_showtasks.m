function varargout = tom_showtasks( jmName, jobName )
%TOM_SHOWTASKS creates visual timing report for distributed tasks
%
%   varargout = plottasks( jmName, jobName )
%
%  PLOTTASKS creates a visual representation of the tasks which comprise a
%  job in the Distributed Computing toolbox. It shows when tasks started
%  and finished and on which worker machine they ran. With no arguments,
%  PLOTTASKS uses the first available job manager it can find and the last
%  finished job in the job manager's list. With one output argument, a
%  handle to the newly created figure is returned. Each worker in the
%  cluster must have a unique name for the display to be correct.
%
%  PLOTTASKS(JOBMANAGER) uses the job manager whose name is JOBMANAGER
%
%  PLOTTASKS(JOBMANAGER,JOB) finds and uses the job whose name is JOB in
%  JOBMANAGER.
%
%EXAMPLE
%   ... = tom_showtasks(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
% Copyright 2005 The MathWorks, Inc.
% Author: Russell.Goyder@mathworks.co.uk
%
%   $Revision: 1.5 $
%   $Date: 2005/05/31 22:02:51 $

% prepare progress bar
wb = waitbar( 0, 'Finding job manager...' );

% grab a job manager
if nargin < 1
    jm = findResource( 'jobmanager' );
else
    jm = findResource( 'jobmanager', 'name', jmName );
end
if isempty( jm )
    delete( wb );
    error('Unable to find a job manager');
end 

% default job manager is the first one we find
jm = jm(1);

% check available
if strcmp( jm.State, 'unavailable')
    delete( wb );
    error('Job manager is unavailable');
else
    waitbar( 0.1, wb, 'Finding job...' );
end 
    
% default job is last job in job manager
if nargin < 2
    jobs = findJob( jm, 'state', 'finished' );
    lj = jobs(end);
else
    lj = findJob( jm, 'name', jobName );
end
if isempty( lj )
    delete( wb );
    error('Job not found');
end

% check it has finished
if ~strcmp( lj.State, 'finished' )
    delete( wb );
    error('Job has not finished');
else
    waitbar( 0.2, wb, 'Getting list of workers' );
end

% get list of workers used for this job
tasks = lj.Tasks;
nTasks = length( tasks );
allWrs = cell( nTasks, 1 );
for i = 1:nTasks,
    waitbar( 0.2 + 0.2*i/nTasks, wb );
    allWrs{i} = tasks(i).Worker.Name;
end
wrs = unique(allWrs);

% hide names of workers for presentation
nWrs = length( wrs );
anonWrs = cell( nWrs, 1 );
for i = 1:nWrs,
   anonWrs{i} = ['Worker ' num2str(i)];
end

% for each task, get start and end times
for i = 1:nTasks,

    waitbar( 0.4 + 0.3*i/nTasks, wb, sprintf( 'Examining task %i', i ) );
    
    st = tasks(i).StartTime;
    ft = tasks(i).FinishTime;
    
    % take off the year info at the end (string cannot be handled directly
    % by datenum?) (and take off day at beginning)
    st(1:4) = [];
    ft(1:4) = [];
    st(end-8:end) = [];
    ft(end-8:end) = [];
    
    % covert to datenum
    tms(i,1) = datenum(st,'mmm dd HH:MM:SS');
    tms(i,2) = datenum(ft,'mmm dd HH:MM:SS');
end
waitbar( 0.7, wb, 'Creating figure' );

% not interested in absolute time, and want seconds
tms = tms - min(tms(:));
tms = tms * 24 * 3600;

% set up figure
f = figure;
ax = handle(axes('parent',f));
ax.YTick = 1:length(wrs);
ax.YLimMode = 'manual';
ax.YLim = [0 length(wrs)+1];
ax.YTickLabel = anonWrs;

% create patches
barHalfWidth = 0.4;
ps = [];
for i = 1:nTasks,

    waitbar( 0.7 + 0.3*i/nTasks, wb, sprintf( 'Creating patch %i', i ) );
    
    % find which worker and generate y coords
    wInd = find( strcmp( wrs, allWrs{i} ) );
    yc = [ wInd-barHalfWidth; wInd-barHalfWidth; ...
        wInd+barHalfWidth; wInd+barHalfWidth ];
    
    % generate x coords based on times
    t1 = tms(i,1);
    t2 = tms(i,2);
    xc = [t1; t2; t2; t1];
    
    ps = [ps; handle(patch(...
        'parent',ax,...
        'xdata', xc,...
        'ydata',yc,...
        'FaceColor',[0.8 0.2 0.2],...
        'EdgeColor','k',...
        'HitTest','off',...
        'visible','on'))];
    
end
waitbar( 1, wb, 'Finishing plot' );

% get serial time
serialTime = sum( tms(:,2) - tms(:,1) );

% get parallel time
parallelTime = max( tms(:) );

% get speedup
speedupFactor = round( serialTime / parallelTime * 10 ) / 10;

% do some labelling
xlabel( 'Time (s)' );
title( [ num2str(nTasks) ' tasks. Serial time: ' num2str( serialTime ) ...
    's. Parallel time: ' num2str( parallelTime ) 's. Speed-up: ' num2str( speedupFactor ) ] );
set( ax.Title, 'interpreter','none' );
set( ax, 'xgrid', 'on' );

% delete waitbar
delete( wb );

% perhaps return handle to the figure
if nargout > 0
    varargout{1} = f;
end
