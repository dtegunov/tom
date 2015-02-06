function [paraell_param]=tom_build_paraell_param(paraell_param);
%TOM_BUILD_PARAELL_PARAM creates ...
%
%   [paraell_param]=tom_build_paraell_param(paraell_param)
%
%PARAMETERS
%
%  INPUT
%   paraell_param       ...
%  
%  OUTPUT
%   paraell_param		...
%
%EXAMPLE
%   ... = tom_build_paraell_param(...);
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

if (isempty(paraell_param)==1)
    paraell_param.jobmanager='default_jobmanager';
    paraell_param.packageloss=0.05;
    paraell_param.number_of_tasks=1;
    paraell_param.workers.min=1;
    paraell_param.workers.max=8;
    paraell_param.timeout=18000;
    paraell_param.restart_workers=0;
end;

if (isfield(paraell_param,'jobmanager')==0)
    paraell_param.jobmanager='default_jobmanager'; 
end;

if (isfield(paraell_param,'packageloss')==0)
     paraell_param.packageloss=0.05;
end;

if (isfield(paraell_param,'workers')==0)
    paraell_param.workers.min=1;
    paraell_param.workers.max=8;
end;

if (isfield(paraell_param.workers,'min')==0)
    paraell_param.workers.min=1;
end;

if (isfield(paraell_param.workers,'max')==0)
    paraell_param.workers.max=1;
end;

if (isstruct('timeout')==0)
    paraell_param.timeout=18000;
end;

if (isstruct('restart_workers')==0)
     paraell_param.restart_workers=0;
end;



