function out=tom_rescale3d(in,newsize,method,waitbarflag)
%TOM_RESCALE3D resizes 3d-data using 2d slices and imresize
%
%   out=tom_rescale3d(in,newsize,method,waitbarflag)
%
%PARAMETERS
%
%  INPUT
%   in                  input volume
%   newsize             new size as a vector
%   method              'nearest'  nearest neighbor interpolation
%                       'bilinear' bilinear interpolation
%                       'bicubic'  bicubic interpolation (default)
%   waitbarflag         flag for showing the waitbar. 0 false, 1 true
%  
%  OUTPUT
%   out                 rescaled volume
%
%EXAMPLE
%   in=tom_emread('Proteasome.vol');
%   out=tom_rescale3d(in.Value,[64 64 64]);
%
%REFERENCES
%
%SEE ALSO
%   imresize
%
%   created by SN 05/30/05
%   updated by AK 10/07/05 waitbar added 
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


if nargin<3
    method='bicubic';
    waitbarflag = 0;
end;

if nargin<4
    waitbarflag = 0;
end

out_tmp=zeros([newsize(1:2) size(in,3)],class(in));

if waitbarflag == 1
    h = waitbar(0,'Rescaling');
    nriterations = size(in,3)+size(out_tmp,1);
    counter = 0;
end


for iz=1:size(in,3)
    if waitbarflag == 1
        waitbar(counter./nriterations, h, ['Rescaling (' num2str(round((counter./nriterations).*100)) ' % done)']);
        counter = counter + 1;
    end
    out_tmp(:,:,iz)=imresize(in(:,:,iz),newsize(1:2),method);
end

out=zeros(newsize,class(in));


for ix=1:size(out_tmp,1)
    if waitbarflag == 1
        waitbar(counter./nriterations, h, ['Rescaling (' num2str(round((counter./nriterations).*100)) ' % done)']);
        counter = counter + 1;
    end
    out(ix,:,:)=imresize(squeeze(out_tmp(ix,:,:)),newsize(2:3),method);
end

if waitbarflag == 1
    close(h);
end
