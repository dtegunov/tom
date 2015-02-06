function zoomed = tom_zoom(vol,ibin)
% TOM_ZOOM zooms
%
%   zoomed = tom_zoom(vol,ibin)
%
%   TOM_ZOOM performs 'inverted binning', i.e. image is magnified WITHOUT
%   interpolation. 
%
% PARAMETERS
%  INPUT
%   VOL         volume
%   IBIN        binning to be inverted
%
%  OUTPUT
%   ZOOMED      2^ibin times magnified volume
%
% EXAMPLE
%
%   mask = tom_spheremask(ones(64,64,64),17,2);
%   zoomed = tom_zoom(mask,1);
%   % -> voxel mask(1,1,1) is assigned to voxels(1:2,1:2,1:2) of zoomed. 
%   tom_dev(tom_bin(zoomed,1)-mask) % difference ~0 !
%
% last change 03/31/05 FF - docu updated
[nx,ny,nz] = size(vol);
nxz = (2^ibin)*nx;nyz = (2^ibin)*ny;nzz = (2^ibin)*nz;
zoomed = zeros(nxz,nyz,nzz);
indx = ndgrid(1:2:nx,1:ny,1:nz);
zoomed(1:2:nxz,1:2:nyz,1:2:nzz) = vol;
zoomed(2:2:nxz,1:2:nyz,1:2:nzz) = vol;
zoomed(1:2:nxz,2:2:nyz,1:2:nzz) = vol;
zoomed(1:2:nxz,1:2:nyz,2:2:nzz) = vol;
zoomed(2:2:nxz,2:2:nyz,1:2:nzz) = vol;
zoomed(1:2:nxz,2:2:nyz,2:2:nzz) = vol;
zoomed(2:2:nxz,1:2:nyz,2:2:nzz) = vol;
zoomed(2:2:nxz,2:2:nyz,2:2:nzz) = vol;