function im=tom_circle(radius);
%TOM_CIRCLE generates a circle
%
%   im=tom_circle(radius);
%
%PARAMETERS
%
%  INPUT
%   filelabel           ...
%   threshold           ...
%   label               ...
%   color               ...
%   transformmatrix     ...
%   iconposition        ...
%   host                ...
%  
%  OUTPUT
%   data		...
%
%   IM=TOM_CIRCLE(radius) generates a circle.It is required only one input, which 
%   is the radius of the circle.The Output of this function is a matrix that contains
%   the generated circle and its dimensions are (2*radius+1)x(2*radius+1).
%
%EXAMPLE
%   im=tom_circle(3)
%
%
%                           0     0     0     1     0     0     0
%                           0     1     1     1     1     1     0
%                           0     1     1     1     1     1     0
%                   im =    1     1     1     1     1     1     1
%                           0     1     1     1     1     1     0
%                           0     1     1     1     1     1     0
%                           0     0     0     1     0     0     0
%
%REFERENCES
%
%SEE ALSO
%   TOM_SPHERE, TOM_ERROR
%
%   created by AL 08/09/02
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

     
center=radius+1;
lims=[1 2*radius+1];
[x,y]=meshgrid(lims(1):lims(2));
im=sqrt((x-center).^2+(y-center).^2)<=radius;