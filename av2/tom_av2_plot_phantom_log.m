function num=tom_av2_plot_phantom_log(log)
%TOM_AV2_PLOT_PHANTOM_LOG creates ...
%
%   num=tom_av2_plot_phantom_log(log)
%
%PARAMETERS
%
%  INPUT
%   log                 ...
%  
%  OUTPUT
%   num                 ...
%
%EXAMPLE
%   ... = tom_av2_plot_phantom_log(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by .. mm/dd/yy
%   updated by ..
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

num=zeros(72,1);

for i=1:size(log,1)
    
    ind=((log(i,2)+180)./5)+1;
    num(ind)=num(ind)+1;    
    
end;


figure; plot(num);