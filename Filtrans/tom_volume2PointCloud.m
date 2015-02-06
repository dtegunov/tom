function cloud=tom_volume2PointCloud(vol,thr,samp)
%tom_volume2PointCloud transform a volume in a pointcloud (set of coordinates)
%
%    cloud=tom_volume2PointCloud(vol,thr,samp)
%
%PARAMETERS
%
%  INPUT
%   vol                    path of the volumes
%   thr                    threshold
%   samp                   sampling of binarized vol 4 final pointcloud
%
%  OUTPUT
%   cloud                  matrix of x y z coordinates
%
%
%EXAMPLE
%
% cloud=tom_volume2PointCloud(tom_spheremask(ones(64,64,64)),0.9,5);
% figure; plot3(cloud(:,1),cloud(:,2),cloud(:,3),'r+');
%
%
%REFERENCES
%
%SEE ALSO
%   
%  tom_align_pointcloud
%
%   created by FB 01/24/12
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


mask=vol>thr;

cloud=zeros((size(vol,1)*size(vol,2)*size(vol,3))./(samp.^3),3);
zz=0;
for ix=1:samp:size(vol,1)
    for iy=1:samp:size(vol,2)
        for iz=1:samp:size(vol,3)
            if (mask(ix,iy,iz) > 0 )
                zz=zz+1;
                cloud(zz,:)=[ix,iy,iz];
            end;
        end;
    end;
end;
cloud=cloud(1:zz,:);


