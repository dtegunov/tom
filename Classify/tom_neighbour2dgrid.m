function [nodes,direction] = tom_neighbour2dgrid(node,gridsize,grid_top)
%TOM_NEIGHBOUR2DGRID creates ...
%
%   [nodes,direction] = tom_neighbour2dgrid(node,gridsize,grid_top)
%
%PARAMETERS
%
%  INPUT
%   node                ...
%   gridsize            ...
%   grid_top            ...
%  
%  OUTPUT
%   nodes        		...
%   direction     		...
%
%EXAMPLE
%   ... = tom_neighbour2dgrid(...);
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

error(nargchk(0, 7, nargin, 'struct'))

if nargin < 3
    metric = 'rect';
end

if node > gridsize(1).*gridsize(2)
    error('No such node.');
end

nodes = [];
direction = [];


if isequal(grid_top,'rect')
    if mod(node,gridsize(1)) > 0
        nodes = [nodes node + 1];
        direction = [direction;-1 0];
    end

    if mod(node,gridsize(1)) ~= 1
        nodes = [nodes node - 1];
        direction = [direction;1 0];
    end

    if node < gridsize(1).*gridsize(2)-gridsize(1)
        nodes = [nodes node + gridsize(1)];
        direction = [direction; 0 1];
    end

    if node > gridsize(1)
        nodes = [nodes node - gridsize(1)];
        direction = [direction; 0 -1];
    end
end

if (isempty(direction))
    error('wrong grid topology');
end;

%direction=direction.*-1;