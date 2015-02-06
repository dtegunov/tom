function vol_out=tom_av2_fourier_plot(Align,iter_num,size_vol,sz_plane)
%TOM_AV2_FOURIER_PLOT creates ...
%
%   vol_out=tom_av2_fourier_plot(Align,iter_num,size_vol,sz_plane)
%
%PARAMETERS
%
%  INPUT
%   Align               ...
%   iter_num            ...
%   size_vol            ...
%   sz_plane            ...
%  
%  OUTPUT
%   vol_out             ...
%
%EXAMPLE
%   ... = tom_av2_fourier_plot(...);
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

if nargin<2
    size_vol=[128 128 128];
end;
    

if nargin<3
    sz_plane=2;
end;

vol_out=zeros(size_vol,'single');

Align=tom_av3_align_sum(Align);

sz=size(Align,2);

start_plane=round((size_vol(1))./2) + 1 -round(sz_plane./2);
stop_plane=round((size_vol(1))./2)+1+round(sz_plane./2);

plane=zeros(size_vol,'single');
plane(:,:,start_plane:stop_plane)=1;

ang=zeros(sz,3);

for i=1:sz 
    ang(i,1)=Align(iter_num,i).Angle.Phi; 
    ang(i,2)=Align(iter_num,i).Angle.Psi; 
    ang(i,3)=Align(iter_num,i).Angle.Theta; 
end;


for i=1:sz 
    tmp=tom_rotate(plane,ang(i,:)); 
    vol_out=vol_out+tmp; 
end;






