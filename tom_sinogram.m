function [ sinogram ] = tom_sinogram( im, anglesteps )

sinogram = zeros(size(im, 1), anglesteps);

for i=0:anglesteps-1
    rotated = imrotate(im, i * 360 / anglesteps, 'bicubic', 'crop');
    rotated = tom_rotate(im, i * 360 / anglesteps);
    rotated = tom_spheremask(rotated, size(im,1)/2);
    proj = sum(rotated)';
    proj = tom_norm(proj,'mean0+1std');
    sinogram(:,i+1) = proj;
end;

end

