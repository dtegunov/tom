
function volume = tom_reconstruction3d(projection_file,myext, marker_file, irefmark, offset,pre_binning,post_binning ...
                    ,weighting,reconstruction_file,filter,reconstruction_size,demo_mode,thickness,r,ireftilt ...
                    ,weightedfile_name,weighted_proj_only,parallel) 
%TOM_RECONSTRUCTION3 performs a weighted 3D-backprojection.
%
% vol = tom_reconstruction3d(projection_file,myext, marker_file, irefmark, offset,pre_binning,post_binning, ...
%       weighting,reconstruction_file,filter,reconstruction_size,demo_mode,thickness,r,ireftilt,weightedfile_name,weighted_proj_only)  
%
%  This function needs a markerfile and the projections of the tiltseries.
%  The projections are moved and rotated before they are weighted. 
%  After the weighting the are backprojected.   
%
%PARAMETERS
%
%   INPUT
%   filelabel                   ...
%    projection_file            name of tiltseries [string]
%    myext                      extension (of file), e.g. '.em' [string]
%    marker_file                name of markerfile [string]
%    irefmark                   index of reference marker 
%    offset[3]                  shift vector of reconstructed volume
%    pre binning                binning of projections before processing (...to speed up the process)
%    post binning               binning of original projections after weighting the projections
%    weighting                  weighting scheme 'none', 'analytical', 'exact'
%    reconstruction_file        name of 3D-reconstruction [string]
%    filter                     radius of low-pass filter in %/100
%    reconstruction_size[3]     dimensions of reconstructed 3D-volume
%    demo_mode                  switches the Visualisation on/off - set to
%                                1/0.
%    thickness                  thickness of reconstructed volume
%                                (required for exact weighting, in pixels)  
%    r[3]                       coordinates assigned to reference marker                                
%    ireftilt                   number of reference tiltangle - default:
%                                0 degree projection
%    weightedfile_name          Name of the weighted file
%    weighted_proj_only         1=creation of weighted file only, 0= weighted file + volume 
%    parallel                   1=parallel mode, 0=simple mode
%
%   OUTPUT
%    volume                     Volume (single!)
%   
%
%   A file ("reconstruction_file") is stored. In the Comment 
%   irefmark, ireftilt, r, offset, binning, filter, weighting
%   are stored.
%
%EXAMPLES
%    (this examples only work if the path tom/data/ is set)
%    vol = tom_reconstruction3d('pyrodictium_','.em', 'mark_pyrodictium.em', 1, [0 0 0],0,1,'exact','testvol.em',1,[512 512 128],0,128,'BMP_');
%    vol = tom_reconstruction3d('pyrodictium_','.em', 'mark_pyrodictium.em', 1, [0 0 0],1,0,'analytical','testvol.em',1,[512 512 128],1,'BMP_');
%    vol = tom_reconstruction3d('pyrodictium_','.em', 'mark_pyrodictium.em', 1, [0 0 0],0,1,'none','testvol.em',1,[512 512 128],0,'BMP_');
%
%REFERENCES
%
%SEE ALSO
%   tom_weight3d,tom_backproj3d,tom_dist
%
%   created by FF 04/17/03
%   updated by FF 04/14/04 analytical weighting tested
%   updated by FB 04/21/04 speed up
%   updated by FF 06/29/04 fixed bug for offset reconstruction
%   updated by WDN 01/04/05 set the name to the weighted file, offset
%   reconstruction
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
error(nargchk(1,18,nargin));
if (nargin < 4)
    irefmark = 1;
end;
if parallel==1
    jobmanager='default_jobmanager';
    jobname='test';
    file=tom_emread(marker_file);
    Matrixmark=file.Value;
    message=tom_alignweight_parallel(jobmanager,jobname, projection_file, myext, marker_file, irefmark, pre_binning, post_binning ...
                    ,weighting,filter,thickness,r,ireftilt,weightedfile_name,size(Matrixmark,2));
    disp(['First message = ' message]);
    message2=tom_backproj3d_parallel(jobmanager, jobname, reconstruction_file, reconstruction_size, size(Matrixmark,2), weightedfile_name, 8);
    disp(['Second message = ' message2]);
    volume=tom_emreadc(reconstruction_file);
else
file=strcat(projection_file,  num2str(1), myext);
file=tom_reademheader(file);
imdim = file.Header.Size(1);
file=tom_emread(marker_file);
Matrixmark=file.Value;
ntilt=size(Matrixmark,2);
hh=tom_emreadc(strcat(projection_file,'1', myext));
sx=size(hh.Value,1);
filter=round(filter*(sx./2));
if ~isempty(findstr(projection_file,'/'))
    a=findstr(projection_file,'/');
    b=size(a);
    data_path=projection_file(1:a(b(2)));
else
    data_path=pwd;
end

if (nargin < 17)
   weighted_proj_only=0;
end;

if isempty(weightedfile_name)
    weightedfile_name='TEMP_BPP_mat_';
end
% determine lowest tiltangle projection
if (nargin < 15)
    [dummy, imintilt]=min(abs(Matrixmark(1,:,1)));
else
    imintilt=ireftilt;
end;
if (nargin < 14)
    %  define r(1), r(2), and r(3) of irefmark
    r = [Matrixmark(2,imintilt,irefmark) Matrixmark(3,imintilt,irefmark) (imdim/2 +1)];
end;
if (nargin < 12)
    demo_mode=0;
end;

% calculate tilt axis and shifts of projection images
[Matrixmark, beta, sigma, x, y, z]  = tom_alignment3d(Matrixmark, irefmark, imintilt, r, imdim);
%betaindeg = beta.*180.0./pi; % perform rotations by beta
betaindeg = beta.*180/pi; % perform rotations by beta
psiindeg = betaindeg+90;  % angle between tilt axis and y-axis
psi=psiindeg*pi/180.0; %tilt axis - determined by alignment program
theta_ang = Matrixmark(1,:,1)'; %projection angle, i.e. read out from the microscope goniometer
% create masks (needed for weighting) outside the loop to speed up the process
[w_func,c_mask]=create_masks(strcat(projection_file,'1', myext),filter,pre_binning);
% calculate offset in not aligned image to display square in demomode
%demo_offset=calc_demo_offset(offset,psi,Matrixmark,imintilt);



 if (pre_binning >=1)
     w_func=tom_bin(w_func,pre_binning);    
     c_mask=tom_bin(c_mask,pre_binning);         
 end;

if (demo_mode ==1)
    demo1;
end;


if (weighted_proj_only==0)
    % allocate some memory
    volume=zeros(reconstruction_size,'single');
end;


for itilt = 1:ntilt,
  
    % read the projection files
    file=strcat(projection_file, num2str(itilt), myext);
    image=tom_emreadc(file);
    orig_header=image.Header;
    if (pre_binning >=1)
        image=tom_bin(image.Value,pre_binning);    
        image=double(image);    
    else
        image=double(image.Value); 
    end;
    
    % normalize to contrast
    imdev = mean(mean(image));  
    image = (image - imdev)./imdev;  %   subtract mean and norm to mean
    
    
    if (demo_mode==1)
        %demo2(image,demo_offset,reconstruction_size,post_binning); %modif.will
        demo2(image,reconstruction_size,pre_binning,post_binning);
    end;
    
    % smooth borders to prevent high contrast oszillations
    image = tom_smooth(image, 12);
    image=tom_rotate(image,-psiindeg,'linear',[floor(size(image,1)./2) floor(size(image,1)./2)]); 
    
    % calculate shifts 
    shift_tx = Matrixmark(7,itilt,1)./(2^pre_binning);
    shift_ty = Matrixmark(8,itilt,1)./(2^pre_binning); 
    
    shiftx = cos(psi)*shift_tx + sin(psi)*shift_ty;
    shifty = -sin(psi)*shift_tx + cos(psi)*shift_ty;
    
    % apply int shifts in real(move) and the rest in fourier space to
    % prevent periodic artifacts
    image=tom_move(image,[-fix(shiftx) -fix(shifty)]);
    fimage = fftshift(fft2(image));
    fimage = shift_in_fs(fimage,[-(shiftx-fix(shiftx)) -(shifty-fix(shifty))] );
    
    % perform weighting
    if isequal(weighting,'exact') | isequal(weighting,'analytical')
        if (isequal(weighting,'analytical'))
            fimage=tom_weight3d('analytical',fimage,filter,'w_func',w_func,'c_mask',c_mask);
        else
            fimage=tom_weight3d('exact',fimage,filter,psiindeg,itilt,theta_ang,thickness); 
        end;
    end;
    image = real(ifft2(ifftshift(fimage))); 
    
    % do the post binning
    if (post_binning >= 1)
        image=tom_bin(image,post_binning);
    end;
    
    %if (weighted_proj_only==0)
        if ( (size(image,1) < reconstruction_size(1)) | (size(image,2) < reconstruction_size(2) ) ) 
            disp('error: reconstruction size bigger than image size');     
            return ;     
        end;
        if (demo_mode==1)
            demo3(image,offset,reconstruction_size);
        end;
    %end;
    
    angle_the = Matrixmark(1,itilt,1);
    angle_phi = 0.0;
    % write out weighted proj
    tmp = strcat(weightedfile_name,num2str(itilt),'.em');
    imt = tom_emheader(image);
    imt.Header=orig_header;
    imt.Header.Size(1)=size(imt.Value,1);
    imt.Header.Size(2)=size(imt.Value,2);
    imt.Header.Size(3)=1; % projection images are two-dimensional
    imt.Header.Objectpixelsize=orig_header.Objectpixelsize.*2.^(pre_binning+post_binning);
    imt.Header.Tiltaxis=0;
    imt.Header.Magic(4)=5; % write floats
    org_path=pwd;
    cd(data_path);
    tom_emwrite(tmp,imt);
    cd(org_path);
    image = imt.Value;%put in to fill header correctly, FF
    
    if (weighted_proj_only==0) %reconstruction
        sx = (size(image,1)-reconstruction_size(1))/2+1;
        sy = (size(image,2)-reconstruction_size(2))/2+1;
        if (offset == [0 0 0])
            image=image(sx:sx+reconstruction_size(1)-1,sy:sy+reconstruction_size(2)-1);
        end;%FF: if set to fix bug for offset ~=0
        image=single(image);
        % do the backprojection       
        off_set=offset./((2^pre_binning)*(2^post_binning));
        tom_backproj3d(volume,image, angle_phi, angle_the, off_set);
        if (demo_mode==1)
            demo4(image,volume,weighting);
        end;
    end;
    disp([' processed micrograph ' num2str(itilt) ' of ' num2str(ntilt)]);
    %when 'Stop' button is clicked in tom_rec3d, stop the procedure
    if ~isempty(findobj(0,'Tag','rec3d'))
        if strcmp('DO reconstruction',get(findobj(0,'Tag','Do_reconstruction'),'String'))%stop the procedure
            %msgbox('Reconstruction stopped','Stop','help');
            break
        end
        
    end       
end;

if (weighted_proj_only==0)

    volume=tom_emheader(volume);
    
    %volume.Value=double(volume.Value);
    volume.Header=orig_header;
    volume.Header.Comment=[num2str(irefmark) ' ' num2str(ireftilt) ' ' num2str(round(r)) ' ' num2str(offset) ' ' num2str(pre_binning) ' ' num2str(post_binning)  ' ' num2str(filter) ' ' weighting];
    volume.Header.Comment=char(volume.Header.Comment');
    
    volume.Header.Size=size(volume.Value);
    volume.Header.Objectpixelsize=orig_header.Objectpixelsize.*2.^(pre_binning+post_binning);
    volume.Header.Magic(4)=5; %write floats
    volume.Header.Tiltaxis=betaindeg;
    volume.Header.Tiltangle=0;
    tom_emwrite(reconstruction_file,volume);
    volume=volume.Value;
else
    volume=-1;
end;
end
%*********************************************
%************  Other Function  ***************
%*********************************************
%----------------------------------------
%------------  create_masks  ------------
%----------------------------------------
function [w_func,c_mask]= create_masks(name,filter,pre_binning)
pic_h=tom_emreadc(name);
dim_h=size(pic_h.Value);
[xh,yh] = ndgrid(-dim_h(1)./2:((dim_h(2)./2)-1));
w_func = tom_norm(abs(xh),1);
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
m=findobj('Tag','r3d');
set(m,'Units','Pixel');pr3d=get(m,'Position');set(m,'Units','Normalized')
if isempty(findobj(0,'Tag','demo1'))
    h_demo=figure('Name','Reconstruction','Tag','demo1'); pdemo=get(h_demo,'Position');
    set(h_demo,'Position', [pr3d(1)+pr3d(3)+10 pr3d(2) pdemo(3) pr3d(4)-50]);
    set(h_demo,'doublebuffer','on','Units','Normalized');
end
drawnow;
%----------------------------------------
%----------------  Demo2  ---------------
%----------------------------------------
function demo2(image,reconstruction_size,pre_binning,post_binning)
h_rec3d=findobj(0,'Tag','rec3d');handles=guidata(h_rec3d);
set(findobj(0,'Tag','Do_reconstruction'),'Enable','off');
m=findobj('Tag','demo1');set(0,'CurrentFigure',m);
ax1=subplot(2,2,1);
set(gca,'position',[0.15 0.581098 0.327023 0.373902]);
tom_imagesc(image,'noinfo','range',[-0.5 0.5]);
set(findobj(0,'Tag','Do_reconstruction'),'Enable','on');
title('IMAGE','FontWeight','bold','Units','normalized');
hold on;
offset(1)=str2num(get(handles.offsetx,'String'))./(2^pre_binning);
offset(2)=str2num(get(handles.offsety,'String'))./(2^pre_binning);;
middle=size(image)./2;
%as(1)=offset(1)./(2^pre_binning);as(2)=offset(2)./(2^pre_binning);
pos=middle-((reconstruction_size(1)./2)*(2^post_binning))+offset;
w=reconstruction_size(1)*(2^post_binning);
h=reconstruction_size(2)*(2^post_binning);
if pos(1)==0 & pos(2)==0
    h_rec=rectangle('Position', [pos(1)+1 pos(2)+1 w h]);
else
    h_rec=rectangle('Position', [pos(1) pos(2) w h]);
end

set(h_rec,'Edgecolor',[1 0 0],...
    'LineWidth',2.5,...
    'Tag','rectan1');
%************ ORIGINAL FLORIAN ***************
%if (post_binning==0 )
%    middle=size(image)./2;
%    pos=middle+offset-(reconstruction_size(1)./2);
%    h=rectangle('Position', [pos(1) pos(2) reconstruction_size(1) reconstruction_size(1)]);
%    set(h,'Edgecolor',[1 0 0]);
%    set(h,'LineWidth',2.5);
%end;
drawnow;
hold off;
%----------------------------------------
%----------------  Demo3  ---------------
%----------------------------------------
function demo3(image,offset,reconstruction_size)
set(findobj(0,'Tag','Do_reconstruction'),'Enable','off');
m=findobj('Tag','demo1');set(0,'CurrentFigure',m);
subplot(2,2,2);
a=round(size(image,1)./2)-(reconstruction_size(1)./2);
if (a==0) a=1; end;
b=round(size(image,1)./2)+(reconstruction_size(1)./2);
if (offset(1)==0 & offset(2)==0)
    m=findobj('Tag','demo1');set(0,'CurrentFigure',m);
    tom_imagesc(image(a:b,a:b),'noinfo');
else
    m=findobj('Tag','demo1');set(0,'CurrentFigure',m);
    tom_imagesc(image,'noinfo');
end;
set(findobj(0,'Tag','Do_reconstruction'),'Enable','on');
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
a=round(size(volume,1)./2);%a=round(size(volume,1)./2)-5;
b=round(size(volume,1)./2);%b=round(size(volume,1)./2)+5;
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

function demo_off=calc_demo_offset(offset,psi,Matrixmark,ref)
demo_off(1)=offset(1);
demo_off(2)=offset(2);

%calculate rotation Matrix
M=[cos(psi) sin(psi) ; -sin(psi) cos(psi)];

% calculate shifts 
shift_tx = Matrixmark(7,ref,1);
shift_ty = Matrixmark(8,ref,1); 
shiftx = cos(psi)*shift_tx + sin(psi)*shift_ty;
shifty = -sin(psi)*shift_tx + cos(psi)*shift_ty;

% transform offset
demo_off=demo_off*M;
demo_off=demo_off-[shiftx shifty];


















