function av3_mergevolume(volfilename, subvolfilename, append, overlap, nx, ny, xdim, ydim, xydim)
% AV3_MERGEVOLUME merges subvolume to large (CCF) volume
%
%   av3_mergevolume(volfilename, subvolfilename, append, overlap, nx, ny, xdim, ydim, xydim)
%   
%PARAMETERS
%   volfilename      : name of resulting volume (em-file)
%   subvolfilename   : name of subvolumes - are assumed as
%                       <subvolfilename>_number.em
%   append           : appendix of subvolumes (e.g. '.em')
%   overlap          : overlap  - diameter of mask (even!)
%   nx               : number of subvol in x
%   ny               : number of subvol in y
%   xdim             : x-dimension of resulting volume
%   ydim             : y-dimension of resulting volume
%   xydim            : dimension in xy (default: 256)
%
% SEE ALSO
%   AV3_CUTVOLUME, oscar
%
%   FF

error(nargchk(1,9,nargin))
if nargin < 9
    xydim = 256;
end;
volume = zeros(xdim, ydim, 256); 
overlap2 = overlap / 2;
x = overlap2 + 1
y = overlap2 + 1;
z = 1;
ind = 1;
for ix=1:nx
    y = 1;
    for iy=1:ny
        filename = [subvolfilename '_' num2str(ind) append];
        subvolume = tom_emread(filename);
        subvolume = tom_red(subvolume.Value, [overlap2+1, overlap2+1, 1], [xydim-overlap, xydim-overlap, 256] );
        volume = tom_paste(volume, subvolume,[x, y, 1]);
        y = y + xydim - overlap;
        ind = ind + 1
    end;
    x = x + xydim - overlap;
end;
tom_emwrite(volfilename, volume);
