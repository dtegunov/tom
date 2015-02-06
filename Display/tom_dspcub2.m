function out=tom_gallery(in,dims,spacings,border,background)
%TOM_GALLERY creates an mxn gallery with spacing and border
%
%   out=tom_gallery(in,dims,spacings,border,background)
%
%PARAMETERS
%
%  INPUT
%   in                  image stack
%   dims                [m n], images in m and n direction
%   spacings            [1 1], spacing in pixels in m and n direction
%   border              [2 2], border in x and y direction
%   background          ...
%  
%  OUTPUT
%   out         		image
%
%EXAMPLE
%   out=tom_gallery(in.Value,[5 8],[3 3],[5 5],-.35);
% 
%   or:
%
%   out=tom_gallery(in.Value,[5 8],[3 3],[5 5],-.35);
%
%REFERENCES
%
%SEE ALSO
%   tom_dspcub, montage
%
%   created by FB, SN date
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

if nargin<5
    background=mean2(in);
end;
if nargin<4
    border=[5 5]; % 5 pixels border
end;
if nargin<3
    spacings=[3 3]; % 3 pixels spacings
end;
if nargin<2
    dims(1)=ceil(sqrt(size(in,1))); 
    dims(2)=ceil(sqrt(size(in,2))); 
end;

m=dims(1);
n=dims(2);
out=ones(size(in,1).*n+(n-1).*spacings(2)+border(2).*2,size(in,2).*m+(m-1).*spacings(1)+border(1).*2).*background;

dim_m=size(in,1);
dim_n=size(in,2);

start_m=border(1)+1;
start_n=border(2)+1;
i_image=0;

for i_m=1:m
    
    for i_n=1:n
        
        i_image=i_image+1;
        
        try
        out(start_n:start_n+dim_n-1,start_m:start_m+dim_m-1)=in(:,:,i_image); % faster
        catch
        end;
        start_n=start_n+dim_n+spacings(2);
        
        
    end
    start_n=border(2)+1;
    start_m=start_m+dim_m+spacings(1);

end

