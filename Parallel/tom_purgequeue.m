function tom_purgequeue(jobmanagername)
%TOM_PURGEQUEUE purges a jobmanager queue
%
%   tom_purgequeue(jobmanagername)
%
%PARAMETERS
%
%  INPUT
%   jobmanagername      name of the jobmanager
%  
%  OUTPUT
%
%EXAMPLE
%   tom_purgequeue(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by AK 01/25/06
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

try
    if nargin == 1
        jm = findResource('scheduler','type','jobmanager','lookupurl',jobmanagername);
    else
        jm = findResource('scheduler','type','jobmanager');
    end

catch
    err = lasterror;
    if strcmp(err.identifier,'distcomp:findResource:LicenseUnavailable') == 1
        error('All distributed computing licenses are in use, please try again later.');
        return;
    else
        error('Unknown error');
        return;
    end
end
finished_jobs = findJob(jm);
if (isempty(finished_jobs)==1)
   disp('Queue is already purged.'); 
   return;
end;
cancel(finished_jobs);
destroy(finished_jobs);

disp('Queue purge completed.');