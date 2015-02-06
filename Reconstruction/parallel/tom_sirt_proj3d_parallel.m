function errorstring=tom_sirt_proj3d_parallel(jobmanager, jobname, name_tiltseries,ext,projections,nr_iterations,binning,size_of_volume,volume_name,diffprefix,nrall_proj)
%TOM_SIRT_PROJ3D_PARALLEL creates ...
%
%   errorstring=tom_sirt_proj3d_parallel(jobmanager, jobname, name_tiltseries,ext,projections,nr_iterations,binning,size_of_volume,volume_name,diffprefix,nrall_proj)
%
%PARAMETERS
%
%  INPUT
%   jobmanager          ...
%   jobname             ...
%   name_tiltseries     ...
%   ext                 ...
%   projections         ...
%   nr_iterations       ...
%   binning             ...
%   size_of_volume      ...
%   volume_name         ...
%   diffprefix          ...
%   nrall_proj          ...
%  
%  OUTPUT
%   errorstring         ...
%
%EXAMPLE
%   ... = tom_sirt_proj3d_parallel(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by ... (author date)
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

tom_emwritec(volume_name, size_of_volume,'new');
unix(['chmod 666 ' volume_name]);
errorstring = '';

for lauf=1:nr_iterations
    disp(['Iteration ' num2str(lauf) ' of ' num2str(nr_iterations)]);
    j = createJob(jm,'Name', jobname);
    set(j,'FileDependencies',{'/fs/bmsan/apps/tom_dev/Parallel/tom_sirt_proj3d_parallel.m'})
    
    for projnr=1:projections
        createTask(j, @tom_sirt_proj3d_worker, 1, {name_tiltseries,ext,projnr,projnr,binning,size_of_volume,volume_name,diffprefix,nrall_proj});
    end

    submit(j);
    waitForState(j);
    out = getAllOutputArguments(j);
    
    if isempty(out)
        tmp = get(j);
        tmp = tmp.Tasks(1);
        errorstring = strcat(errorstring, 'on worker "', tmp.Worker.Hostname, '":  ', tmp.Errormessage, ', all workers failed with this error.');
        destroy(j);
        return;
    end;
    
    for i=1:size(out,1)
        if isempty(out{i})
            tmp = get(j);
            tmp = tmp.Tasks(i);
            errorstring = strcat(errorstring, 'on worker "', tmp.Worker.Hostname, '":  ', tmp.Errormessage);
            destroy(j);
            return;
        end;
    end;
    destroy(j);
    errorstring=[errorstring tom_backproj3d_parallel(jobmanager, jobname, volume_name, size_of_volume, projections, diffprefix, 8)];
end;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function status=tom_sirt_proj3d_worker(name_tiltseries,ext,min_nr_proj,max_nr_proj,binning,size_of_volume,volume_name,diffprefix,nrall_proj)

addpath('/fs/bmsan/apps/tom_dev/Reconstruction/');
addpath('/fs/bmsan/apps/tom_dev/IOfun/');
addpath('/fs/bmsan/apps/tom_dev/Parallel/');
addpath('/fs/bmsan/apps/tom_dev/Sptrans/');

volume=tom_emreadc(volume_name);
sx=size(volume,1);
sx_2=size(volume,1)./2;

%mask=tom_spheremask(ones(size(volume)),sx_2,1,[sx_2+1 sx_2+1 sx_2+1]);
%mask_proj=tom_spheremask(ones(sx,sx),sx_2,3,[sx_2+1 sx_2+1 1]);

for i=min_nr_proj:max_nr_proj
    proj=tom_emreadc([name_tiltseries num2str(i) ext]);
    proj.Value=tom_bin(single(proj.Value),binning);
    proj.Value=(proj.Value-mean2(proj.Value))./mean2(proj.Value);
    new_proj=tom_proj3d(single(volume.Value),[0 proj.Header.Tiltangle]);
    new_proj=new_proj./(250.*nrall_proj.*(1./cos(proj.Header.Tiltangle.*pi./180)).*sx);
    tmp=single(proj.Value)-single(new_proj);
    tmp=tom_emheader(tmp);
    tmp.Header.Tiltangle=proj.Header.Tiltangle;
%    tmp.Header.Size(1)=proj.Header.Size(1)./(2.^binning);
%    tmp.Header.Size(2)=proj.Header.Size(2)./(2.^binning);
%    tom_emwrite([diffprefix num2str(i) '.em'], tmp);
    tom_emwrite([diffprefix num2str(i) '.em'], tmp);
    %tom_emwrite('/tmp/parallel.tmp',tmp);
    %unix(['mv /tmp/parallel.tmp ' diffprefix num2str(i) '.em']);
end;
status = 1;


