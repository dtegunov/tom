function mtf = tom_createmtf(image)
%TOM_FILTER convolutes with spherical or quadratic kernel
%
%  mtf = tom_createmtf(image)
%
%PARAMETERS
%
%  INPUT
%   image               edge image (edge has to be upright)
%  
%  OUTPUT
%   mtf                 resulting mtf (1D)
%
%EXAMPLE
%
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created 01/09/07 AK and AZ
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

edgepos = zeros(1,size(image,1));

%for i=301:2048-100
for i=1:size(image,1)
   line = image(:,i);
   line = smooth(line,10);

   linea = line(1:end-6);
   lineb = line(7:end);
   linediff = linea - lineb;

   pos = find(line(4:end-3)-linediff< 0);
   edgepos(i) = pos(1)+3;
   
end

imnew = zeros(size(image,1)./2,size(image,2));
%imnew = zeros(257,2048-500);
%for i=301:2048-100
    
for i=1:size(image,1)
    line = image(:,i);
    tmp = edgepos(i)-size(image,1)./4+1:edgepos(i)+size(image,1)./4;
%    tmp = tmp(find(tmp>0));
    line = line(tmp);
    %line = line(edgepos(i)-128:edgepos(i)+128);
    imnew(:,i)=line;
end

imnew=tom_taper(imnew,[2048 2048]);
avgedge = mean(imnew,2);

psf = zeros(1,length(avgedge)-1);
for i=2:length(avgedge)
   psf(i-1) = avgedge(i-1)-avgedge(i); 
end

mtf = abs(fft(psf));
mtf = mtf(1:size(image,1)./4);
%mtf = mtf(1:256);
mtf = smooth(mtf,10);
mtf = tom_norm(mtf,1);

