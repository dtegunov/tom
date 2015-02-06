function all_pos=tom_codevecs2distmap(codvecs, gridsize, grid_top,metric,obs_dim,origin)
%TOM_CODEVECS2DISTMAP creates ...
%
%   all_pos=tom_codevecs2distmap(codvecs, gridsize, grid_top,metric,obs_dim,origin)
%
%PARAMETERS
%
%  INPUT
%   codvecs             ...
%   gridsize            ...
%   grid_top            ...
%   metric              ...
%   obs_dim             ...
%   origin              ...
%  
%  OUTPUT
%   all_pos     		...
%
%EXAMPLE
%   .. = tom_codevecs2distmap(...);
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

error(nargchk(0, 6, nargin, 'struct'))

if (nargin < 3)
    grid_top='rect';
end;


if (nargin < 4)
    metric='euclidian';
end;


if (nargin < 5)
    obs_dim=1;
end;


if (nargin < 6)
    origin=[(floor(gridsize(1)./2)+1) (floor(gridsize(1)./2)+1)];
end;


%normalize all values to phase contrast
for i=1:size(codvecs,1)
    codvecs(i,:)=tom_norm(codvecs(i,:)+100,'phase');
end;

%build nodes structure
ns = struct();
for i=1:gridsize(1).*gridsize(2)
    dists = [];
    [neighbours directions] = tom_neighbour2dgrid(i,gridsize,grid_top);
    ns(i).neighbours = neighbours;
    ns(i).directions = directions;
    ns(i).pos='';
    for ii=1:length(neighbours)
        if (strcmp(metric,'euclidian') )
            dists(ii) = sum( ( (codvecs(i,:)-codvecs(neighbours(ii),:) ).^2)) .^(1/2);
        else
            %reshape images
            if (size(obs_dim,2)==2)
                im_node=reshape(codvecs(i,:),obs_dim(1),obs_dim(2));
                im_neighbour=reshape(codvecs(neighbours(ii),:),obs_dim(1),obs_dim(2));
            end;

            if (size(obs_dim,2)==3)
                im_node=reshape(codvecs(i,:),obs_dim(1),obs_dim(2),obs_dim(3));
                im_neighbour=reshape(codvecs(neighbours(ii),:),obs_dim(1),obs_dim(2),obs_dim(3));
            end;

            ccf=tom_corr(im_node,im_neighbour,'norm');
            [pos val]=tom_peak(ccf);
            dists(ii)=-0.5.*val + 0.5;
        end;
    end;
    ns(i).dists = dists;
end


%initialize
nodepos=((origin(1)-1).*gridsize(1)) + origin(2);
ns(nodepos).pos =[0 0];


while all_pos_set(ns)

    %find nodes with positions
    z=1;
    for i=1:size(ns,2)
        if (isempty(ns(i).pos)==0)
            list(z)=i;
            z=z+1;
        end;
    end;

    %find neighbours & calc distance
    for i=1:max(size(list))
        tmp_neig=ns(list(i)).neighbours;
        for ii=1:max((size(tmp_neig)))
            if (isempty(ns(tmp_neig(ii)).pos))
                ns=calc_pos(ns,tmp_neig(ii));
            end;
        end;
    end;



end;

%build output strucut
for i=1:size(codvecs,1)
    all_pos(i,1) = ns(i).pos(1);
    all_pos(i,2) = ns(i).pos(2);
end;



%helper functions
function ns=calc_pos(ns,node)

tmp_neig=ns(node).neighbours;

new_pos=[0 0];

for i=1:max(size(tmp_neig))
    if (isempty(ns(tmp_neig(i)).pos )==0 )
        new_pos=new_pos+ (ns(node).directions(i,:).*ns(node).dists(i)) + (abs(ns(node).directions(i,:)).*ns(tmp_neig(i)).pos);
    end;
end;

%update
ns(node).pos=new_pos;


function [flag]=all_pos_set(ns)


flag=0;

for i=1:size(ns,2)
    if (isempty(ns(i).pos))
        flag=1;
    end;
end


