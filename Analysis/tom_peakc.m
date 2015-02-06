function [c val in] =tom_peakc(in,radius)
%TOM_tom_PEAKC creates ...
%
%   [c val in] =tom_peakc(in,radius)
%
%PARAMETERS
%
%  INPUT
%   in                  ...
%   radius              ...
%  
%  OUTPUT
%   c           		...
%   val           		...
%   in           		...
%
%EXAMPLE
%   ... tom_amira_createisosurface(...);
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

if nargin==1
    radius=1;
end;

if (size(in,1)==1)
    s=1;
    dims(1)=int32(size(in,2));
else
    s=size(size(in),2);
    dims=int32(size(in));
end;

cut_size=[2.*radius 2.*radius 2.*radius];

if (s==2)
    dims(3)=1;
    cut_size=[2.*radius 2.*radius];
end

if (s==1)
    dims(3)=1;
    dims(2)=1;
    cut_size=[2.*radius];
end;

c=int32([0 0 0]);
val=single(0);
in=single(in);

tom_peakinc(in,dims,val,c);
out=tom_cut_out(in,(double(c)-double(radius)),cut_size,'nofill');
mask= ones(size(out));
mask = tom_spheremask(mask,radius,0,[size(out,1)./2+1 size(out,2)./2+1 size(out,3)./2+1]);

 mask=(mask==0);
out=double(out).*mask;
if (s>1)
    in=tom_paste(in,out,(double(c)-double(radius)));
    
else
    start_p=double(c(1)-(size(mask,2)./2));
    stop_p=double(start_p+size(mask,2))-1;
    in(start_p:stop_p)=0;
end;
c=double(c);
val=double(val);




