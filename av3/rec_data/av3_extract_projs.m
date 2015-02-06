function projstack = av3_extract_projs(motl,tomofilename,direction,dims,nlayers)
% AV3_EXTRACT_PROJS projects subtomograms
%
%   projstack = av3_extract_projs(motl,tomofilename,direction,dimensions,nlayers)
%
%   The purpose of AV3_EXTRACT_PROJS is the generation of projections of
%   the aligned tomograms to have an idea about their SNR etc.
%
% PARAMETERS
%  INPUT
%   MOTL         motivelist
%   TOMOFILENAME name of tomogram (EM format)
%   direction    direction of projection (1=x, 2=y, 3=z)
%   dimensions   dimensions of projection (square!)
%   nlayers      number of layers to be projected (around center)
%
%  OUTPUT
%   PROJSTACK    stack containing the different projections (3rd coordinate
%                   = particleindex)
%
% last change FF 03/31/05 - updated doc

%tomo = tom_emread(tomofilename);
if direction ==1
    projstack = zeros(dims(2),dims(3),size(motl,2));
elseif direction ==2
    projstack = zeros(dims(1),dims(3),size(motl,2));
elseif irection == 3
    projstack = zeros(dims(1),dims(2),size(motl,2));
else
    disp('choose direction = 1, 2 or 3 !');
    exit;
end;
for ii=1:size(motl,2)
    x= motl(8,ii);y= motl(9,ii);z= motl(10,ii);
    x=x-floor(dims(1)/2);y=y-floor(dims(2)/2);z=z-floor(dims(3)/2);
    part = tom_emreadc(tomofilename, 'subregion', [x y z],[dims(1)-1 dims(2)-1 dims(3)-1]);
    phi=motl(17,ii);
    psi=motl(18,ii);
    theta=motl(19,ii);
    part = tom_spheremask(double(tom_rotate(part.Value,[-psi,-phi,-theta])),floor(dims(1)/2)-2,1);
    if direction ==1
        proj = sum(part(floor((dims(1)-nlayers)/2)+1:floor((dims(1)-nlayers)/2)+nlayers,:,:));
    elseif direction ==2
        proj = sum(part(:,floor((dims(2)-nlayers)/2)+1:floor((dims(2)-nlayers)/2)+nlayers,:));
    elseif irection == 3
        proj = sum(part(:,:,floor((dims(3)-nlayers)/2)+1:floor((dims(3)-nlayers)/2)+nlayers));
    else
        disp('choose direction = 1, 2 or 3 !');
        exit;
    end;
    projstack(:,:,ii) = squeeze(proj);
    %part = tom_red(vol.Value,[x-32,y-32,z-32],[64,64,64]);
end;
