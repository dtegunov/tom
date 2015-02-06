function tom_show_cluster_status_local(flag)
%TOM_SHOW_CLUSTER_STATUS start gkrellm system monitors
%
%   tom_show_cluster_status(flag)
%
%PARAMETERS
%
%  INPUT
%   flag                ('off') turns all gkrellm off
%  
%  OUTPUT
%
%EXAMPLE
%   tom_show_cluster_status;
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


if nargin==0
    flag='on';
end;
list={'pasing','sendling','hadern','haidhausen','halifax','calgary','solln','schwabing','giesing'};
if strcmp(flag,'on')
    x_coord=10;
    for i=1:length(list)
        unix(['ssh -Y ' list{i} ' gkrellm -g +' num2str(x_coord) '+25 &']); 
        x_coord=x_coord+240;
        pause(1);
    end;
    
else
    for i=1:length(list)
        unix(['ssh ' list{i} ' killall gkrellm &']); 
        pause(.1);
    end;
end;
