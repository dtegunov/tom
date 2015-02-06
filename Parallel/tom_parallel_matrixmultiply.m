function c = tom_parallel_matrixmultiply(a,b,parallelstruct)
%TOM_PARALLEL_MATRIXMULTIPLY creates ...
%
%   c = tom_parallel_matrixmultiply(a,b,parallelstruct)
%
%PARAMETERS
%
%  INPUT
%   a                   ...
%   b                   ...
%   parallelstruct      ...
%  
%  OUTPUT
%   c                   ...
%
%EXAMPLE
%   ... = tom_parallel_matrixmultiply(...);
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

if size(a,2) ~= size(b,1)
    error('The number of columns of a must equal the number of rows of b');
end

if nargin < 3
    parallelstruct = tom_parallelsettings();
end

size_c = [size(a,1),size(b,2)];

%write input matrices
tom_emwritec('/fs/bmsan/pub/mtimes_a.em',a,'standard','single');
tom_emwritec('/fs/bmsan/pub/mtimes_b.em',b,'standard','single');

%create output matrix
tom_emwritec('/fs/bmsan/pub/mtimes_c.em',[size_c 1],'new','single');

unix('chmod 777 /fs/bmsan/pub/mtimes_a.em /fs/bmsan/pub/mtimes_b.em /fs/bmsan/pub/mtimes_c.em');

numtasks = parallelstruct.number_of_tasks;
jm = findResource('jobmanager','name',parallelstruct.jobmanager);
j = createJob(jm,'Name','autopicker');
pdd={'/fs/bmsan/apps/tom_dev/Parallel/tom_parallel_matrixmultiply.m' '/fs/bmsan/apps/tom_dev/IOfun/'};
set(j,'FileDependencies',pdd);
packages=tom_calc_packages(numtasks,size(a,2));

disp('Sending jobs to workers...');
numtasks2 = numtasks;
for i=1:numtasks
    range = [packages(i,1),packages(i,2),packages(i,3)];
    createTask(j,@tom_parallel_matrixmultiply_worker,0,{size_c,range});
end
disp('Processing on workers...');
submit(j);
waitForState(j);
out = getAllOutputArguments(j);
[result_p errorsum]=tom_get_paraell_results(j);
if (errorsum > 0);
    for i=1:size(result_p,1)
        disp(['Hostname: ' result_p{i,1} '  Error: ' result_p{i,2}]);
    end;
end;
destroy(j);
disp('Finished');
c = tom_emreadc('/fs/bmsan/pub/mtimes_c.em');
c = c.Value;

delete('/fs/bmsan/pub/mtimes_a.em');
delete('/fs/bmsan/pub/mtimes_b.em');
delete('/fs/bmsan/pub/mtimes_c.em');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tom_parallel_matrixmultiply_worker(size_c,range)

a = tom_emreadc('/fs/bmsan/pub/mtimes_a.em');
h = tom_reademheader('/fs/bmsan/pub/mtimes_b.em');
b = tom_emreadc('/fs/bmsan/pub/mtimes_b.em','subregion',[1 range(1) 1],[h.Header.Size(1)-1 range(3)-1 0]);
a = a.Value;
b = b.Value';

c = zeros([size_c(2),range(2)-range(1)+1],'single');

for i=1:size(c,1)
    for j=1:size(c,2)
        c(i,j) = sum(a(i,:).*b(j,:));
    end
end

tom_emwritec('/fs/bmsan/pub/mtimes_c.em',c,'subregion',[1 range(1) 1],[size(c,1) size(c,2) 1]);