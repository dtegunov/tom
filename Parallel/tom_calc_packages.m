function packages=tom_calc_packages(num_of_nodes,num_of_calculations,index)
%TOM_CALC_PACKAGES creates ...
%
%   packages=tom_calc_packages(num_of_nodes,num_of_calculations,index)
%
%PARAMETERS
%
%  INPUT
%   num_of_nodes        number of cpus
%   num_of_calculations number of calculations
%   index               ...
%  
%  OUTPUT
%   packages            ...
%
%EXAMPLE
%   ... = tom_calc_packages(...);
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



if (nargin==2)
    index=ones(num_of_nodes,1);
end;

package_size=floor(num_of_calculations./num_of_nodes);
package_rest=mod(num_of_calculations,num_of_nodes);

start=1;
zz=1;
for i=1:num_of_nodes
    if ((num_of_nodes-i+1)==package_rest)
        package_size=package_size+1;
    end;
    if (index(i)==1)
        packages(zz,:)=[start (start+package_size-1) package_size];
        zz=zz+1;
    end;
    start=start+package_size;
end;

