function tom_display_som(classes,inputvect,im_size,codevect)
%TOM_DISPLAY_SOM creates ...
%
%   tom_display_som(classes,inputvect,im_size,codevect)
%
%PARAMETERS
%
%  INPUT
%   classes             ...
%   inputvect           ...
%   im_size             ...
%   codevect            ...
%  
%  OUTPUT
%
%EXAMPLE
%   tom_display_som(...);
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

error(nargchk(0, 4, nargin, 'struct'))




%build codevect stack
for i=1:size(codevect,1)
    codev_tmp(:,:,i)=reshape(codevect(i,:),im_size);
end;
figure; tom_dspcub(codev_tmp); set(gcf,'Name','Code Vectors');

avg_tmp=zeros(im_size(1),im_size(2),size(codevect,1));
avg_num=zeros(size(codevect,1),1);

%build class averages
for i=1:size(inputvect,1)
    avg_tmp(:,:,classes(i))=avg_tmp(:,:,classes(i))+reshape(inputvect(i,:),im_size);
    avg_num(classes(i))=avg_num(classes(i))+1;
end;

%normalize it don't critisize it !
for i=1:size(codevect,1)
     avg_tmp(:,:,i)=avg_tmp(:,:,i)./avg_num(i);
end;

figure; tom_dspcub(avg_tmp); set(gcf,'Name','Class Averages');
