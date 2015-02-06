function tom_cov_paralell(infilename,outfilename,paraell)
%TOM_COV_PARALLEL creates ...
%
%   tom_cov_paralell(infilename,outfilename,paraell)
%
%PARAMETERS
%
%  INPUT
%   infilename          ...
%   outfilename         ...
%   paraell             ...
%  
%  OUTPUT
%
%EXAMPLE
%   tom_cov_paralell(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by ... (author date)
%   updated by ...
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

error(nargchk(0, 3, nargin, 'struct'))


h=tom_reademheader(infilename);
sz_m=h.Header.Size;
num_of_tasks=paraell.number_of_tasks;
%max_lost_packages=paraell.packageloss;
jobmanager=paraell.jobmanager;



tom_emwritec('out.em',[sz_m(2) sz_m(2) 1],'new');

[a,b]=system(['chmod -R ugo+rwx ' infilename]);
[a,b]=system(['chmod -R ugo+rwx ' outfilename]);

jm = findResource('jobmanager','name',jobmanager);
% set(jm.jobs,'MaximumNumberOfWorkers',40);
result_p=tom_get_paraell_results(jm,'hosts');
disp(['classifying on ' num2str(size(result_p,1)) ' hosts']);
for iii=1:size(result_p,1)
    disp(result_p{iii});
end;
clear('result_p');
j = createJob(jm,'Name','demojob');

pdd={'/fs/bmsan/apps/tom_dev/IOfun' '/fs/bmsan/apps/tom_dev/Classify'};
%     '/fs/bmsan/apps/tom_dev/Filtrans/'  '/fs/bmsan/apps/tom_dev/Analysis/' '/fs/bmsan/apps/tom_dev/Misc/' '/fs/bmsan/apps/tom_dev/Sptrans/' ...
%     '/fs/bmsan/apps/tom_dev/Geom/' '/fs/bmsan/apps/tom_dev/av2/'  };
set(j,'FileDependencies',pdd)
pdd2={[infilename] [outfilename]};
set(j,'PathDependencies',pdd2);
rehash toolbox;

%set(j,'FileDependencies',{'/fs/bmsan/apps/tom_dev/av2/tom_av2_multi_ref_alignment.m'})

packages=tom_calc_packages(num_of_tasks,sz_m(2));
for i=1:num_of_tasks
    createTask(j,@tom_cov_paralell_worker,1,{infilename,outfilename,packages(i,1:2),i});
end;
submit(j);
%tom_disp_paraell_progress(j,packages(:,3));
waitForState(j);
out = getAllOutputArguments(j);
[result_p errorsum]=tom_get_paraell_results(j);
if (errorsum > 0);
    for i=1:size(result_p,1)
        disp(['Hostname: ' result_p{i,1} '  Error: ' result_p{i,2}]);
    end;
end;
destroy(j);

disp('end');






