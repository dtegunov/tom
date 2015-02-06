function zoomed = tom_zoom(vol,ibin)
%
%   zoomed = tom_zoom(vol,ibin)
%
%   TOM_ZOOM performs 'inverted binning', i.e. image is magnified WITHOUT
%   interpolation. Example:
%   mask = tom_spheremask(ones(64,64,64),17,2);
%   zoomed = tom_zoom(mask,1);
%   % -> voxel mask(1,1,1) is assigned to voxels(1:2,1:2,1:2) of zoomed. 
%   tom_dev(tom_bin(zoomed,1)-mask) % difference ~0 !
%
% 17/01/05 modified by GS
%   - works now for 2d and 3d
%   - ibin > 1 allowed
%   
%
% 

[nx,ny,nz] = size(vol);
switch nz
    case 1
        for ind =1:ibin
            nxz = (2^ind)*nx;nyz = (2^ind)*ny;
            zoomed = zeros(nxz,nyz,1);
            %indx = ndgrid(1:2:nx,1:ny,1:nz);
            zoomed(1:2:nxz,1:2:nyz) = vol;
            zoomed(2:2:nxz,1:2:nyz) = vol;
            zoomed(1:2:nxz,2:2:nyz) = vol;
            zoomed(2:2:nxz,2:2:nyz) = vol;
            vol = zoomed;
            ind = ind+1;
        end
    otherwise
        for ind =1:ibin
            nxz = (2^ind)*nx;nyz = (2^ind)*ny;nzz = (2^ind)*nz;
            zoomed = zeros(nxz,nyz,nzz);
            %indx = ndgrid(1:2:nx,1:ny,1:nz);
            zoomed(1:2:nxz,1:2:nyz,1:2:nzz) = vol;
            zoomed(2:2:nxz,1:2:nyz,1:2:nzz) = vol;
            zoomed(1:2:nxz,2:2:nyz,1:2:nzz) = vol;
            zoomed(1:2:nxz,1:2:nyz,2:2:nzz) = vol;
            zoomed(2:2:nxz,2:2:nyz,1:2:nzz) = vol;
            zoomed(1:2:nxz,2:2:nyz,2:2:nzz) = vol;
            zoomed(2:2:nxz,1:2:nyz,2:2:nzz) = vol;
            zoomed(2:2:nxz,2:2:nyz,2:2:nzz) = vol;
            vol = zoomed;
            ind = ind+1;
        end
end
 