function im = av3_vol2gallery(vol,nx,idspcubmode,mask,flag)
%  AV3_VOL2GALLERY can be used for converting a 3D volume into a 2D gallery 
%
%   im = av3_vol2gallery(vol,nx,idspcubmode,mask,flag)
%
%   vol         volume (3D)
%   nx          number of frames plotted into x-dimension
%   idspcubmode dspcubmode - 0,1,2 for xy, xz, and yz slices
%   mask        mask (optional) - is multiplied before gallery is created,
%               i.e. it has to be 3D with dimensions of VOL
%   flag        flag (optional) - can be set to 'norm' for normalizing VOL
%               to standard deviation.
%
%   im          2D output (as in tom_dspcub)
%
%   VOL2GALLERY converts 3D volume into a 2D gallery. This function is
%   thought to be used as an interface for the EM program to use its
%   (2D)-classification capabilities.
%
%   FF 30/06/04
if nargin <3
    idspcubmode = 0;
end;
if nargin>4
    if isequal(flag,'norm')
        [mv mn mx stv] = tom_dev(vol,'noinfo');
        vol = (vol-mv)/stv;
    end;
end;
if nargin>3
    vol = vol.*mask;
end;
switch idspcubmode
    case 0
        ix = 1;
        iy=1;
        for ii=1:size(vol,3)
            x = (ix-1)*size(vol,1)+1;y = (iy-1)*size(vol,2)+1;
            im(x:x+size(vol,1)-1,y:y+size(vol,2)-1) = vol(:,:,ii);
            ix = ix+1;
            if ix > nx
                ix =1;
                iy = iy+1;
            end;
        end;
    case 1
        ix = 1;
        iy=1;
        for ii=1:size(vol,2)
            x = (ix-1)*size(vol,1)+1;y = (iy-1)*size(vol,3)+1;
            im(x:x+size(vol,1)-1,y:y+size(vol,3)-1) = squeeze(vol(:,ii,:));
            ix = ix+1;
            if ix > nx
                ix =1;
                iy = iy+1;
            end;
        end;
    case 2
        ix = 1;
        iy=1;
        for ii=1:size(vol,1)
            x = (ix-1)*size(vol,1)+1;y = (iy-1)*size(vol,3)+1;
            im(x:x+size(vol,2)-1,y:y+size(vol,3)-1) = squeeze(vol(ii,:,:));
            ix = ix+1;
            if ix > nx
                ix =1;
                iy = iy+1;
            end;
        end;
end;
