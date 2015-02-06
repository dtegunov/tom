function av3_cutvolume(volfilename, subvolfilename, overlap, nx, ny, xydim)
% AV3_CUTVOLUME        cuts large volume into pieces with overlap
%
%   av3_cutvolume(volfilename, subvolfilename, overlap, nx, ny,xydim)
%   
%PARAMETERS
%   volfilename      : name of volume (em-file)
%   subvolfilename   : name of subvolumes - will be stored as
%                       <subvolfilename>_number.em
%   overlap          : overlap  - diameter of mask
%   nx               : number of subvol in x
%   ny               : number of subvol in y
%   xydim            : dimension in xy (default: 256)
%
%   FF
error(nargchk(1,6,nargin))
if nargin < 6
    xydim = 256;
end;
volume = tom_emreadc(volfilename);volume = double(volume.Value);
xdim = size(volume,1); ydim = size(volume,2); zdim = size(volume,3);
%head = tom_reademheader(volfilename);
%xdim = double(head.Header.Size(1));ydim = double(head.Header.Size(2));zdim = double(head.Header.Size(3))

x = 1;
y = 1;
z = 1;
ind = 1;
for ix=1:nx
    y = 1;
    for iy=1:ny
        %xxx = tom_emreadc(volfilename,'subregion',[x, y, z], [xydim-1, xydim-1, zdim-1]);
        %subvolume=zeros(xydim, xydim, 256);
        %subvolume(:,:,1:zdim) = double(xxx.Value);
        subvolume=tom_red(volume,[x, y, z], [xydim, xydim, 256]);
        filename = [subvolfilename '_' num2str(ind) '.em'];
        tom_emwrite(filename, subvolume);
        y = y + xydim - overlap;
        ind = ind + 1
    end;
    x = x + xydim - overlap;
end;
