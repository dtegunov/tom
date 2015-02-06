function errorstring = tom_alignweight_parallel(jobmanager,jobname,projection_file,myext, marker_file, irefmark, pre_binning, post_binning ...
                    ,weighting,filter,thickness,r,ireftilt,weightedfile_name,tilts)
%TOM_ALIGNWEIGHT_PARALLEL performs an alignment and optional weighting.
%
% errorstring = tom_alignweight_parallel(jobmanager, jobname, projection_file,myext, marker_file, irefmark, pre_binning,post_binning, ...
%       weighting,filter,thickness,r,ireftilt,weightedfile_name,tilts)  
%
%  This function needs a markerfile and the projections of the tiltseries.
%  The projections are moved and rotated before they are weighted. 
%
%  INPUT
%   jobmanager              name of jobmanager
%   jobname                 some descriptive name for the job
%   projection_file         name of tiltseries [string]
%   myext                   extension (of file), e.g. '.em' [string]
%   marker_file             name of markerfile [string]
%   irefmark                index of reference marker 
%   pre binning             binning of projections before processing (...to speed up the process)
%   post binning            binning of original projections after weighting the projections
%   weighting               weighting scheme 'none', 'analytical', 'exact'
%   filter                  radius of low-pass filter in %/100
%   reconstruction_size[3]  dimensions of reconstructed 3D-volume
%   thickness               thickness of reconstructed volume
%                            (required for exact weighting, in pixels)  
%   r[3]                    coordinates assigned to reference marker                                
%   ireftilt                number of reference tiltangle - default:
%                            0 degree projection
%   weightedfile_name       Name of the weighted file
%   tilts                   number of projections
%  OUTPUT
%   errorstring             status message (empty means no error, otherwise error is shown)
%   
%
%  All the projections are converted and stored. The working directory
%  MUST be acessible by all the workers
%
%EXAMPLE
%   ... = tom_alignweight_parallel(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   tom_weight3d,tom_backproj3d,tom_dist, tom_rec3d, tom_reconstruction
%
%   created by AK 08/18/05
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
                
jm = findResource('jobmanager','name',jobmanager);
j = createJob(jm,'Name', jobname);
set(j,'FileDependencies',{'/fs/bmsan/apps/tom_dev/Parallel/tom_alignweight_parallel.m'});

for starttilt=1:tilts
    createTask(j, @tom_alignweight_worker, 1, {projection_file,myext, marker_file, irefmark, pre_binning,post_binning ...
                    ,weighting,filter,thickness,r,ireftilt,weightedfile_name,starttilt,1});
end

submit(j);
waitForState(j);

%do some error checking
out = getAllOutputArguments(j);
errorstring = '';
if isempty(out)
    tmp = get(j);
    tmp = tmp.Tasks(1);
    errorstring = strcat(errorstring, 'on worker "', tmp.Worker.Hostname, '":  ', tmp.Errormessage, ', all workers failed with this error.');
    destroy(j);
    return;
end;

for i=1:size(out,1)
    if out{i} ~= 1
        tmp = get(j);
        tmp = tmp.Tasks(i);
        errorstring = strcat(errorstring, 'on worker', tmp.Worker.Hostname, ': ', tmp.Errormessage, '\n');
        destroy(j);
        return;
    end;
end;

destroy(j);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                

function volume = tom_alignweight_worker(projection_file,myext, marker_file, irefmark, pre_binning,post_binning ...
                    ,weighting,filter,thickness,r,ireftilt,weightedfile_name,starttilt, numtilt)
                
addpath('/fs/bmsan/apps/tom_dev/Reconstruction/');
addpath('/fs/bmsan/apps/tom_dev/IOfun/');
addpath('/fs/bmsan/apps/tom_dev/Sptrans/');
addpath('/fs/bmsan/apps/tom_dev/Filtrans/');
addpath('/fs/bmsan/apps/tom_dev/Parallel/');

%unix(['cd ' directory]);
%   check arguments

if (nargin < 4)
    irefmark = 1;
end;
file=tom_reademheader([projection_file '1' myext]);
imdim = file.Header.Size(1);
file=tom_emreadc([marker_file]);
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

if isempty(weightedfile_name)
    weightedfile_name='TEMP_BPP_mat_';
end
% determine lowest tiltangle projection
if isempty(ireftilt)
    [dummy, imintilt]=min(abs(Matrixmark(1,:,1)));
else
    imintilt=ireftilt;
end;
if isempty(r)
    %  define r(1), r(2), and r(3) of irefmark
    r = [Matrixmark(2,imintilt,irefmark) Matrixmark(3,imintilt,irefmark) (imdim/2 +1)];
end;

% calculate tilt axis and shifts of projection images
[Matrixmark, beta, sigma, x, y, z]  = tom_alignment3d(Matrixmark, irefmark, imintilt, r, imdim);
%betaindeg = beta.*180.0./pi; % perform rotations by beta
betaindeg = beta.*180/pi; % perform rotations by beta
psiindeg = betaindeg+90;  % angle between tilt axis and y-axis
psi=psiindeg*pi/180.0; %tilt axis - determined by alignment program
theta_ang = Matrixmark(1,:,1)'; %projection angle, i.e. read out from the microscope goniometer
% create masks (needed for weighting) outside the loop to speed up the process
%cd (directory);
[w_func,c_mask]=create_masks(strcat(projection_file, '1' , myext),filter,pre_binning);




 if (pre_binning >=1)
     w_func=tom_bin(w_func,pre_binning);    
     c_mask=tom_bin(c_mask,pre_binning);         
 end;

for itilt = starttilt:starttilt+numtilt-1,
  
    % read the projection files
    file=strcat(projection_file, num2str(itilt), myext);
    image=tom_emread(file);
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
            fimage=tom_weight3d('exact',fimage,filter,psi,itilt,theta_ang,thickness); 
        end;
    end;
    image = real(ifft2(ifftshift(fimage))); 
    
    % do the post binning
    if (post_binning >= 1)
        image=tom_bin(image,post_binning);
    end;
    
    
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
    %tom_emwrite('/tmp/parallel.tmp',imt);
    %unix(['mv /tmp/parallel.tmp ' tmp]);
    %cd(org_path);
    %image = imt.Value;%put in to fill header correctly, FF
      
           
end;
volume = 1;

%*********************************************
%************  Other Function  ***************
%*********************************************
%----------------------------------------
%------------  create_masks  ------------
%----------------------------------------
function [w_func,c_mask]= create_masks(name,filter,pre_binning)
pic_h=tom_emread(name);
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

