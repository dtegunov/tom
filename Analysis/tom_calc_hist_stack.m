function hist_stack=tom_calc_hist_stack(filename,hist_spacing)
%TOM_CALC_STACK creates ...
%
%   hist_stack=tom_calc_hist_stack(filename,hist_spacing)
%
%PARAMETERS
%
%  INPUT
%   filename            ...
%   hist_spacing        ...
%  
%  OUTPUT
%   hist_stack			...
%
%EXAMPLE
%   ...=tom_calc_hist_stack(...);
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

error(nargchk(0, 2, nargin, 'struct'))

if isstruct(filename)
    align2d = filename;
    for i=1:size(align2d,2)
        particle=tom_emreadc(align2d(1,i).filename,'subregion',[align2d(1,i).position.x-align2d(1,i).radius align2d(1,i).position.y-align2d(1,i).radius 1],[2*align2d(1,i).radius-1 2*align2d(1,i).radius-1 0]);
        hist_stack(:,i)=hist(reshape(particle.Value,1,[]),hist_spacing);
    end;
else

    h=tom_reademheader(filename);
    sz=h.Header.Size;

    for i=1:sz(3)
        tmp=tom_emreadc(filename,'subregion',[1 1 i],[sz(1)-1 sz(2)-1 0]);
        hist_stack(:,i)=hist(reshape(tmp.Value,1,[]),hist_spacing);
    end;
end
hist_stack = hist_stack';