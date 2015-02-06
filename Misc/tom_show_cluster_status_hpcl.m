function tom_show_cluster_status_hpcl(flag,start,stop)
%TOM_SHOW_CLUSTER_STATUS_HPCL start gkrellm system monitors
%
%   tom_show_cluster_status_hpcl(flag,start,stop)
%
%PARAMETERS
%
%  INPUT
%   flag                ('on') turns all gkrellm on, ('off') turns all
%                       gkrellm off
%   start               start # of first HPCL20XX cluster computers.
%   stop                stop # of last HPCL20XX cluster computers.
%
%  maximum of 16 computers can be displayed.
%
%  OUTPUT
%
%EXAMPLE
%   tom_show_cluster_status_hpcl('on',1,16);
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
    start=1;
    stop=16;
end;
if stop-start>16
    stop=start+16;
end;
for i=1:9
    list{i}=['hpcl100' num2str(i)];
end;
for i=10:37
    list{i}=['hpcl10' num2str(i)];
end;
if strcmp(flag,'on')
    x_coord=10;
    for i=start:stop
        unix(['ssh -X -Y ' list{i} ' gkrellm -g +' num2str(x_coord) '+25 --geometry 50 300 &']); 
        x_coord=x_coord+240;
        pause(1);
    end;

end;
if strcmp(flag,'off')
     for i=start:stop
        unix(['ssh ' list{i} ' killall gkrellm &']); 
        pause(.1);
    end;
end;
