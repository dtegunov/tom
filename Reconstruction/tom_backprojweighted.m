
function vol = tom_backprojweighted(projection_file,myext, marker_file, offset,post_binning, ...
                    reconstruction_file,reconstruction_size,filter) 
%TOM_BACKPROJWEIGHTED performs a 3D-backprojection
%
%   vol = tom_backprojweighted(projection_file,myext, marker_file,offset,post_binning, ...
%                    reconstruction_file,reconstruction_size,filter) 
%
%TOM_BACKPROJWEIGHTED performs a 3D-backprojection of previously aligned
%   and weighted projections
%
%  This function needs a markerfile and the projections of the tiltseries.
%  The projections are moved and rotated before they are weighted. 
%  After the weighting the are backprojected.   
%
%PARAMETERS
%
%  INPUT
%   projection_file         name of tiltseries [string] - weighted and
%                            aligned projections
%   myext                   extension (of file), e.g. '.em' [string]
%   marker_file             name of markerfile [string]
%   offset[3]               shift vector of reconstructed volume
%   post binning            binning of original projections after weighting
%                           the projections
%   reconstruction_file     name of 3D-reconstruction [string]
%   reconstruction_size[3]  dimensions of reconstructed 3D-volume
%   filter                  ...
%  
%  OUTPUT
%   vol                     Volume (single!)
%
%   A file ("reconstruction_file") is stored. In the Comment 
%   irefmark, ireftilt, r, offset, binning, filter, weighting
%   are stored.
%
%EXAMPLE
%   (this examples only work if the path tom/data/ is set)
%   vol = tom_reconstruction3d('pyrodictium_','.em', 'mark_pyrodictium.em', 1, [0 0 0],0,1,'none','testvol.em',1,[512 512 128],0);
%
%REFERENCES
%
%SEE ALSO
%   tom_weight3d,tom_backproj3d,tom_dist,tom_recontsruction3d
%
%   created by FF 04/17/03
%   updated by FF 04/14/04 analytical weighting tested
%   updated by FB 04/21/04 speed up
%   updated by FF 06/29/04 fixed bug for offset reconstruction
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
                
%   check arguments
error(nargchk(7,8,nargin));
file=strcat(projection_file,  num2str(1), myext);
file=tom_reademheader(file);
imdim = file.Header.Size(1);
file=tom_emread(marker_file);
Matrixmark=file.Value;
ntilt=size(Matrixmark,2);
hh=tom_emreadc(strcat(projection_file,'1', myext));
sx=size(hh.Value,1);
if (nargin>7)
    filter=round(filter*(sx./2));
else
    filter=(sx./2);
end;
% calculate tilt axis and shifts of projection images
%[Matrixmark, beta, sigma, x, y, z]  = tom_alignment3d(Matrixmark, irefmark, imintilt, r, imdim);
%betaindeg = beta.*180.0./pi; % perform rotations by beta
%psiindeg = betaindeg-90;  % angle between tilt axis and y-axis
%psi=psiindeg*pi/180.0; %tilt axis - determined by alignment program
theta_ang = Matrixmark(1,:,1)'; %projection angle, i.e. read out from the microscope goniometer
% create masks (needed for weighting) outside the loop to speed up the process
% allocate some memory
volume=single(zeros(reconstruction_size));
for itilt = 1:ntilt,
    % read the projection files
    file=strcat(projection_file, num2str(itilt), myext);
    image=tom_emreadc(file);image=image.Value;
    % do the post binning
    if nargin>7
        image = tom_bandpass(image,0,filter);
    end;
    if (post_binning >= 1)
        image=tom_bin(image,post_binning);
    end;
    if ( (size(image,1) < reconstruction_size(1)) | (size(image,2) < reconstruction_size(2) ) ) 
        disp('error: reconstruction size bigger than image size');     
        return ;     
    end;
    angle_the = Matrixmark(1,itilt,1);
    angle_phi = 0.0;
    sx = (size(image,1)-reconstruction_size(1))/2+1;
    sy = (size(image,2)-reconstruction_size(2))/2+1;
    if (offset == [0 0 0])
        image=image(sx:sx+reconstruction_size(1)-1,sy:sy+reconstruction_size(2)-1);
    end;%FF: if set to fix bug for offset ~=0
    image=single(image);
    % do the backprojection
    tom_backproj3d(volume,image, angle_phi, angle_the, offset);
    disp([' processed micrograph ' num2str(itilt) ' of ' num2str(ntilt)]);
end;

volume=tom_emheader(volume);
%here are the most important parameters stored in order of appearance
%volume.Header.Comment=['Parameter: ' num2str(irefmark) ' ' num2str(ireftilt) ' ' num2str(round(r)) ' ' num2str(offset) ' ' num2str(pre_binning) ' ' num2str(post_binning)  ' ' num2str(filter) ' ' weighting]; 
%volume.Value=double(volume.Value);
tom_emwrite(reconstruction_file,volume);
vol=double(volume.Value); 
msgbox('Reconstruction done','Finished','help');
%*********************************************
%************  Other Function  ***************
%*********************************************
%----------------------------------------
%------------  create_masks  ------------
%----------------------------------------
function [w_func,c_mask]= create_masks(name,filter,pre_binning)
pic_h=tom_emreadc(name);
dim_h=size(pic_h.Value);
% if (pre_binning >=1)
%     dim_bin=size(tom_bin(pic_h.Value,pre_binning),1);
% else
%     dim_bin=size(dim_h,1);
% end

% cut_size=(dim_bin)./2;
[xh,yh] = ndgrid(-dim_h(1)./2:((dim_h(2)./2)-1));
w_func = tom_norm(abs(xh),1);

% if (pre_binning >=1)
%     w_func= w_func(  ((dim_h(1)./2+1)-cut_size):((dim_h(1)./2)+cut_size),((dim_h(2)./2+1)-cut_size):((dim_h(2)./2)+cut_size) );
% end
c_mask = sqrt((xh).^2 + (yh).^2) <= filter;

%----------------------------------------
%------------  shift_in_fs   ------------
%----------------------------------------
function im_out=shift_in_fs(im,v)
[dimx,dimy]=size(im);
[x,y]=ndgrid( -floor(size(im,1)/2):-floor(size(im,1)/2)+(size(im,1)-1),...
    -floor(size(im,2)/2):-floor(size(im,2)/2)+size(im,2)-1);
v = v./[dimx dimy];
x = v(1)*x + v(2)*y;clear y; 
im_out=(im.*exp(-2*pi*i*x));

%----------------------------------------
%----------------  Demo1  ---------------
%----------------------------------------
% demo functions for Visualisation
function demo1()
%figure; set(gcf,'Position', [40 50 1200 905]);
m=findobj('Tag','r3d');
set(m,'Units','Pixel');pr3d=get(m,'Position');set(m,'Units','Normalized')
h_demo=figure('Name','Reconstruction','Tag','demo1'); pdemo=get(h_demo,'Position');
set(h_demo,'Position', [pr3d(1)+pr3d(3)+10 pr3d(2) pdemo(3) pr3d(4)-50]);
set(h_demo,'doublebuffer','on','Units','Normalized');
%subplot(2,2,1);
%subplot(2,2,2);
%subplot(2,2,3);
drawnow;

%----------------------------------------
%----------------  Demo2  ---------------
%----------------------------------------
function demo2(image,offset,reconstruction_size,post_binning)
m=findobj('Tag','demo1');set(0,'CurrentFigure',m);
subplot(2,2,1);
m=findobj('Tag','demo1');set(0,'CurrentFigure',m);
set(gca,'position',[0.15 0.581098 0.327023 0.373902]);
tom_imagesc(image,'noinfo','range',[-0.5 0.5]);
title('IMAGE','FontWeight','bold','Units','normalized');
hold on;

if (offset(1)==0 & offset(2)==0 & post_binning==0 )

a=round(size(image,1)./2)-(reconstruction_size(1)./2);
b=round(size(image,1)./2)+(reconstruction_size(2)./2);

h=line([a a],[a b]);
set(h,'color',[1 0 0]);
set(h,'LineWidth',2.5);
drawnow; 
h=line([a b],[b b]);
set(h,'color',[1 0 0]);
set(h,'LineWidth',2.5);
drawnow; 
h=line([b b],[b a]);
set(h,'color',[1 0 0]);
set(h,'LineWidth',2.5);
drawnow; 
h=line([b a],[a a]);
set(h,'color',[1 0 0]);
set(h,'LineWidth',2.5);
drawnow; 
end;
drawnow;    
hold off;    

%----------------------------------------
%----------------  Demo3  ---------------
%----------------------------------------
function demo3(image,offset,reconstruction_size)
m=findobj('Tag','demo1');set(0,'CurrentFigure',m);
subplot(2,2,2);
a=round(size(image,1)./2)-(reconstruction_size(1)./2);
if (a==0) a=1; end;
b=round(size(image,1)./2)+(reconstruction_size(1)./2);
if (offset(1)==0 & offset(2)==0)
    tom_imagesc(image(a:b,a:b),'noinfo');
else
    tom_imagesc(image,'noinfo');
end;
title('WEIGHTED IMAGE','FontWeight','bold','Units','normalized');
m=findobj('Tag','demo1');set(0,'CurrentFigure',m);
set(gca,'position',[0.54 0.581098 0.327023 0.373902]);
drawnow;    


%----------------------------------------
%----------------  Demo4  ---------------
%----------------------------------------
function demo4(image,volume,wei)
m=findobj('Tag','demo1');set(0,'CurrentFigure',m);
subplot(2,2,3);
m=findobj('Tag','demo1');set(0,'CurrentFigure',m);
set(gca,'position',[0.33 0.05 0.327023 0.5]);
axis image;
a=round(size(volume,3)./2)-5;
b=round(size(volume,3)./2)+5;

t_vol=volume(:,:,a:b);
ixy=sum(t_vol,3);
a=round(size(volume,1)./2)-5;
b=round(size(volume,1)./2)+5;
t_vol=volume(a:b,:,:);
iyz=(squeeze(sum(t_vol,1))');   
a=round(size(volume,1)./2)-5;
b=round(size(volume,1)./2)+5;
t_vol=volume(:,a:b,:);
ixz=(squeeze(sum(t_vol,2)));
all=zeros(size(ixy,1)+20+size(iyz,1),size(ixy,2)+20+size(ixz,2));
all(1:size(ixy,1),1:size(ixy,2))=ixy(:,:);
all(size(ixy,1)+21:size(ixy,1)+20+size(iyz,1),1:size(iyz,2))=iyz(:,:);
all(1:size(ixz,1),size(ixy,2)+21:size(ixy,2)+20+size(ixz,2))=ixz(:,:);
if (strcmp(wei,'none')==1)
    tom_imagesc(all,'noinfo','range');  
else
    tom_imagesc(all,'noinfo','range',[-0.25 0.25]);   
end;
title('VOLUME','FontWeight','bold','Units','normalized');
drawnow;



