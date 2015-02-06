function [param cloud_model_alg]=tom_align_pointCloud(cloud_tmpl,cloud_model)
%tom_align_pointCloud aligns 2 pointclouds by rotation and translation 
%
%    [param cloud_model_alg]=tom_align_pointCloud(cloud_tmpl,cloud_model)
%
%PARAMETERS
%
%  INPUT
%   cloud_tmpl                  template cloud
%   cloud_model                 model cloud (will be aligned in respect 2 template)
%   
%
%  OUTPUT
%   param                  alignment parameters shift and rotation matrix
%   cloud_model_alg        aligned point cloud
%
%EXAMPLE
%   
%  test_vol=(tom_spheremask(ones(64,64,64),6,0,[44 44 33])+tom_cylindermask(ones(64,64,64),12))>0;
%  test_vol(:,:,1:15)=0; test_vol(:,:,end-15:end)=0;
%  test_vol_trans=tom_move(tom_rotate(test_vol,[-0.3 0 0]),[0 0 0]);
%  cloud=tom_volume2PointCloud(test_vol,0.87,1);
%  cloud_trans=tom_volume2PointCloud(test_vol_trans,0.87,1);
%  [param cloud_model_alg]=tom_align_pointCloud(cloud,cloud_trans);
%  figure; plot3(cloud(:,1),cloud(:,2),cloud(:,3),'r+'); hold on;  plot3(cloud_trans(:,1),cloud_trans(:,2),cloud_trans(:,3),'go');  
%  figure; plot3(cloud(:,1),cloud(:,2),cloud(:,3),'r+'); hold on;  plot3(cloud_model_alg(:,1),cloud_model_alg(:,2),cloud_model_alg(:,3),'go');
%  disp(' ');
%
%
%REFERENCES
%
%SEE ALSO
%   
%  tom_align_pointCloud
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

if (nargin < 3)
    do_refl=false;
end;

if (nargin < 4)
    do_scale=false;
end;

if (size(cloud_tmpl,1) > size(cloud_model,1))
   cloud_tmpl=cloud_tmpl(1:size(cloud_model,1),:);
else
  cloud_model=cloud_model(1:size(cloud_tmpl,1),:);  
end;

[D cloud_model_alg transform]=procrustes(cloud_tmpl,cloud_model,'Reflection',do_refl,'Scaling',do_scale);

param.transform=transform;
param.D=D;



