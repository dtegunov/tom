function errorstring=tom_backproj3d_parallel(jobmanager, jobname, volname, volsize,numproj, projprefix, numworkers)    
%TOM_ALIGNWEIGHT_PARALLEL performs an alignment and optional weighting
%
%   errorstring=tom_backproj3d_parallel(jobmanager, jobname, volname, volsize,numproj, projprefix, numworkers)
%
%   Performs backprojection in parallel mode like tom_backproj3d. The working directory
%   MUST be acessible by all the workers
%
%PARAMETERS
%
%  INPUT
%   jobmanager          name of jobmanager
%   jobname             some descriptive name for the job
%   volname             full path to the output volume
%   volsize[3]          size of the reconstructed volume
%   numproj             total number of projections
%   projprefix          full path to projections
%   numworkers          number of workers to use (power of 2 is
%                        strongly recommended for speed!)
%  
%  OUTPUT
%   errorstring         status message (empty means no error, otherwise
%                        error is shown)
%
%EXAMPLE
%   ... = tom_backproj3d_parallel(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   tom_weight3d,tom_backproj3d,tom_dist, tom_rec3d, tom_reconstruction,
%   tom_alignweight_parallel
%
%   created by AK 08/18/05
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

jm = findResource('jobmanager','name',jobmanager);
j = createJob(jm,'Name', jobname);
set(j,'FileDependencies',{'/fs/bmsan/apps/tom_dev/Parallel/tom_backproj3d_parallel.m'})

tom_emwritec(volname, volsize,'new');
unix(['chmod 666 ' volname]);

slicez = volsize(3) / numworkers;

for z_offset=0:slicez:volsize(3)-slicez
    createTask(j, @tom_backproj3d_worker, 1, {z_offset, volname, volsize, slicez, numproj, projprefix});
end

submit(j);
waitForState(j);
out = getAllOutputArguments(j);
errorstring = '';
if isempty(out)
    tmp = get(j);
    tmp = tmp.Tasks(1);
    errorstring = strcat(errorstring, 'on worker "', tmp.Worker.Hostname, '":  ', tmp.Errormessage, ', all workers failed with this error.');
    destroy(j);
    return;
end;

for i=1:size(out,1)
    if out{i} ~= 1
        tmp = get(j);
        tmp = tmp.Tasks(i);
        errorstring = strcat(errorstring, 'on worker', tmp.Worker.Hostname, ': ', tmp.Errormessage, '\n');
        destroy(j);
        return;
    end;
end;
destroy(j);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function status=tom_backproj3d_worker(z_offset, volname, volsize, slicez, numproj, projprefix)

addpath('/fs/bmsan/apps/tom_dev/Reconstruction/');
addpath('/fs/bmsan/apps/tom_dev/IOfun/');
addpath('/fs/bmsan/apps/tom_dev/Parallel/');

vol=zeros(volsize(1),volsize(2),slicez,'single'); 

z_offset2=z_offset - volsize(3)/2 + slicez/2;

for n=1:numproj
          proj=tom_emreadc([projprefix num2str(n) '.em']);
          tom_backproj3dc(vol,proj.Value,0,proj.Header.Tiltangle,[0 0 z_offset2]);
end;

tom_emwritec(volname,vol,'subregion',[1 1 z_offset+1],[volsize(1) volsize(2) slicez]);
status = 1;
%status=z_offset;