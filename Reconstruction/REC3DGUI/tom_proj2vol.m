function vol = tom_proj2vol(proj,angles,dimvol,offset)
%%%%%
%%%%%
% proj
% dimvol = [NX NY NZ]
% angles = [phi theta]
% 
% offset = [x y z]

% C implementation
vol = tom_proj2volc(single(proj),double(angles),double(dimvol),double(offset));
  
 
