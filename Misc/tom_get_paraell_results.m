function [result error_sum]=tom_get_paraell_results(j,flag)
%TOM_GET_PARAELL_RESULTS creates ...
%
%   [result error_sum]=tom_get_paraell_results(j,flag)
%
%PARAMETERS
%
%  INPUT
%   j                   ...
%   flag                ...
%  
%  OUTPUT
%   result              ...
%   error_sum           ...
%
%EXAMPLE
%   .. = tom_get_paraell_results(...);
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

if (nargin==1)
    flag='error';
end;
 error_sum=0;

if (strcmp(flag,'error'));

   
    for i=1:size(get(j.Task),1)
        result{i,1}=get(j.Task(i).Worker,'Hostname');
        result{i,2}=get(j.Task(i),'Errormessage');
        if (isempty(get(j.Task(i),'Errormessage')) == 0)
            error_sum=error_sum +1;
        end;
        if (strcmp(result{i,2},'Task cancelled by user'))
            result{i,3}='1';
        else
            result{i,3}='0';
        end;

    end;

else

   num_of_workers=get(j,'NumberOfIdleWorkers');
    
   for i=1:num_of_workers
         result{i,1}=get(j.idleWorkers(i),'Hostname');
   end;
   
end;