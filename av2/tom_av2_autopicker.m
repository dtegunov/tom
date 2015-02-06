function [align2d,ccfmap] = tom_av2_autopicker(template, imagecell, outalign, mask_template, angles_start, angles_stop, angles_incr, number, binning, fftfilename, reusefftfile, maskcell,ccfmapflag,filterstruct,parallelstruct)
%TOM_AV2_AUTOPICKER picks particles using multiple templates on a stack of 2d images
%
%   [align2d, hist_stack] = tom_av2_autopicker(template, imagecell, outalign, mask_template, angles_start, angles_stop, angles_incr, number, binning, parallelstruct)
%
%PARAMETERS
%
%  INPUT
%   template            2d or 3d matrix with template(s)
%   imagecell           cell with complete file and pathnames of images (these
%                        files have to be world readable when using parallel mode)
%   outalign            full pathname and filename to output alignment file
%   mask_template       2d mask matrix (optional, to disable, use [])
%   angles_start        rotational search angle start
%   angles_stop         rotational search angle stop
%   angles_incr         rotational search angle increment
%   number              number of particles to pick from each image
%   binning             apply n times binning before search
%   fftfilename         full path to the file where the ffts of the template
%                       will be stored
%   reusefftfile        use existing fftfile
%   parallelstruct      structure with parallel processing information
%                        (optional) (for more information see
%                        tom_parallelsettings)
%  
%  OUTPUT
%   align2d             2d alignment structure
%   hist_stack          histogram stack of particles for use in tom_pcagui
%
%EXAMPLE
%   tom_amira_createisosurface(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by AK 03/08/06
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


if nargin < 15
    parallelstruct = [];
end

if nargin < 15
    filterstruct = [];
end

if nargin < 13
    ccfmapflag = 0;
end

if nargin < 12
    maskcell = [];
end

if nargin < 9 
    binning = 1;
end

if nargin < 8
    number = 50;
end
if nargin < 7
    angles_incr = 5;
end
if nargin < 6
    angles_stop = 355;
end
if nargin < 5
    angles_start = 0;
end
if nargin < 4
    mask_template = [];
end

disp('Preprocessing...');

if binning > 0
    new_temp = zeros(floor(size(template,1)./2^binning),floor(size(template,2)./2^binning),size(template,3),'single');
    for i=1:size(template,3)
        new_temp(:,:,i) = tom_binc(template(:,:,i),binning);
    end
    template = new_temp;
end

template = single(template);
templatesize = size(template);
angs = [angles_start:angles_incr:angles_stop];

if isempty(mask_template)
    mask_template = tom_spheremask(ones(size(template,1),size(template,2)),size(template,1)./2-2);
end

mask_template = single(mask_template);
[temp mea st P]=norm_inside_mask(template(:,:,1),mask_template);

%generate ffts of template
image = tom_reademheader(imagecell{1});
imagesize = [image.Header.Size(1) image.Header.Size(2)];

if reusefftfile == 0 || reusefftfile == 2
    
    %prepare ffts of template
    filesize = (((imagesize(1)./2^binning).* (imagesize(2)./2^binning) .* (length(angs).*size(template,3))) .*4 + 512)./1024./1024./1024;
    disp(['file size of FFT files will be ' sprintf('%0.3g',2.*filesize) ' GB']);
    disp('creating FFT files, please be patient...')

    tom_emwritec([fftfilename '_real'],[imagesize(1)./2^binning imagesize(2)./2^binning length(angs).*size(template,3)],'new','single');
    tom_emwritec([fftfilename '_img'],[imagesize(1)./2^binning imagesize(2)./2^binning length(angs).*size(template,3)],'new','single');
    unix(['chmod 777 ' fftfilename '_real']);
    unix(['chmod 777 ' fftfilename '_img']);    
    
    if ~isempty(parallelstruct) && ~isempty(fieldnames(parallelstruct))
        numtasks = parallelstruct.number_of_tasks;
        jm = findResource('scheduler','type','jobmanager','lookupurl',parallelstruct.jobmanager);
        j = createJob(jm,'Name','autopicker');
        pdd={'/fs/pool/pool-bmsan-apps/tom_dev/Sptrans' '/fs/pool/pool-bmsan-apps/tom_dev/Analysis/' '/fs/pool/pool-bmsan-apps/tom_dev/IOfun/' '/fs/pool/pool-bmsan-apps/tom_dev/av2/' '/fs/pool/pool-bmsan-apps/tom_dev/Filtrans/' '/fs/pool/pool-bmsan-apps/tom_dev/Geom/'};
        set(j,'FileDependencies',pdd);
        packages=tom_calc_packages(numtasks,size(template,3));
        startidx = 1;
        for i=1:numtasks
            templatestack = template(:,:,packages(i,1):packages(i,2));
            if ~isempty(templatestack)
                createTask(j,@tom_av2_autopicker_createffts_worker,0,{templatestack,mask_template,angs,imagesize,fftfilename,binning,startidx});
                startidx = startidx + size(templatestack,3);
            else
                numtasks = numtasks-1;
            end
        end 
        disp('Processing FFT of template on workers...');
        submit(j);
        waitForState(j);
        [result_p errorsum]=tom_get_paraell_results(j);
        if (errorsum > 0);
            for i=1:size(result_p,1)
                disp(['Hostname: ' result_p{i,1} '  Error: ' result_p{i,2}]);
            end;
            return;
        end;
        destroy(j);
        
    else
        disp('Creating ffts of template on localhost...');
        tom_av2_autopicker_createffts_worker(template,mask_template,angs,imagesize,fftfilename,binning,1);
    end

    if reusefftfile == 2
        align2d = '';
        return;
    end
    
end


if ~isempty(parallelstruct) && ~isempty(fieldnames(parallelstruct))
    numtasks = parallelstruct.number_of_tasks;
    jm = findResource('scheduler','type','jobmanager','lookupurl',parallelstruct.jobmanager);
    j = createJob(jm,'Name','autopicker');
    pdd={'/fs/pool/pool-bmsan-apps/tom_dev/Sptrans' '/fs/pool/pool-bmsan-apps/tom_dev/Analysis/' '/fs/pool/pool-bmsan-apps/tom_dev/IOfun/' '/fs/pool/pool-bmsan-apps/tom_dev/av2/' '/fs/pool/pool-bmsan-apps/tom_dev/Filtrans/' '/fs/pool/pool-bmsan-apps/tom_dev/Geom/'};
    set(j,'FileDependencies',pdd);
    packages=tom_calc_packages(numtasks,size(imagecell,2));
end



%parallel execution
if ~isempty(parallelstruct) && ~isempty(fieldnames(parallelstruct)) && length(imagecell) > 1
    disp('Sending pick jobs to workers...');
    numtasks2 = numtasks;
    for i=1:numtasks
        imagestack = imagecell(packages(i,1):packages(i,2));
        if ~isempty(imagestack)
            createTask(j,@tom_av2_autopicker_worker,1,{imagestack,templatesize,mask_template,angs,number,binning,fftfilename,P,maskcell,0,filterstruct});
        else
            numtasks2 = numtasks2-1;
        end
    end
    disp('Processing on workers...');
    submit(j);
    waitForState(j);
    out = getAllOutputArguments(j);
    [result_p errorsum]=tom_get_paraell_results(j);
    if (errorsum > 0);
        for i=1:size(result_p,1)
            disp(['Hostname: ' result_p{i,1} '  Error: ' result_p{i,2}]);
        end;
    end;

    disp('Processing results...');
    
    if numtasks > 0
    lauf=1;
    for i=1:numtasks
        for k=1:size(out{i},2)
            align2d(lauf) = out{i}(k);
            lauf=lauf+1;
        end
    end
    end
    destroy(j);

%local execution
else
    disp('Processing on localhost...');
    [align2d,ccfmap] = tom_av2_autopicker_worker(imagecell,templatesize,mask_template,angs,number,binning,fftfilename,P,maskcell,ccfmapflag,filterstruct);
end

if ~isempty(outalign)
    save(outalign, 'align2d');
    disp(['Alignment file written to ' outalign]);
end

% disp('Generating histogram stack...');
% 
% hist_bins = 50;
% hist_stack = zeros(hist_bins,size(align2d,2));
% for i=1:size(align2d,2)
%     particle=tom_emreadc(align2d(1,i).filename,'subregion',[align2d(1,i).position.x-align2d(1,i).radius align2d(1,i).position.y-align2d(1,i).radius 1],[2*align2d(1,i).radius-1 2*align2d(1,i).radius-1 0]);
%     hist_stack(:,i)=hist(reshape(double(particle.Value),[],1),hist_bins);
% end;
% tom_emwrite('hist_stack.em',hist_stack');
% 
% disp('Histogram stack written to hist_stack.em. FINISHED.')

disp('FINISHED.');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tom_av2_autopicker_createffts_worker(template,mask_template,angs,imagesize,fftfilename,binning,startidx)

sx_ref_vol_2=floor(imagesize(1)./(2.*2^binning));sx_ref_2=floor(size(template,1)./2);
sy_ref_vol_2=floor(imagesize(2)./(2.*2^binning));sy_ref_2=floor(size(template,2)./2);

lauf = startidx;
for j=1:size(template,3)
    for ang = angs
        temp = tom_rotate(template(:,:,j),ang,'linear');
        temp = single(temp .* mask_template);
        [temp mea st P]=norm_inside_mask(temp,mask_template);
        temp_volume=(zeros(imagesize(1)./(2^binning),imagesize(2)./(2^binning),'single'));
        temp_volume(sx_ref_vol_2-sx_ref_2+1:sx_ref_vol_2+sx_ref_2,sy_ref_vol_2-sy_ref_2+1:sy_ref_vol_2+sy_ref_2)=single(temp);
        tom_emwritec([fftfilename '_real'],single(real(conj(fft2(temp_volume)))),'subregion',[1 1 lauf],[imagesize(1)./2^binning imagesize(2)./2^binning 1]);
        tom_emwritec([fftfilename '_img'],single(imag(conj(fft2(temp_volume)))),'subregion',[1 1 lauf],[imagesize(1)./2^binning imagesize(2)./2^binning 1]);
        lauf = lauf + 1;
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [align,ccfmap] = tom_av2_autopicker_worker(imagestack,templatesize,mask_template,angs,number,binning,fftfilename,P,maskcell,ccfmapflag,filterstruct)

if length(templatesize) == 2
    templatesize(3) = 1;
end

se = strel('disk',20);

align = struct();

image = tom_reademheader(imagestack{1});

sx_ref_vol_2=floor(image.Header.Size(1)./(2.*2^binning));sx_ref_2=floor(templatesize(1)./2);
sy_ref_vol_2=floor(image.Header.Size(2)./(2.*2^binning));sy_ref_2=floor(templatesize(2)./2);


disp('Processing images');
%loop over images
nofiles = length(imagestack);
l = 1;
for filename = imagestack
    
    image = tom_emreadc(filename{1},'binning',binning);
    image = single(image.Value);

    %filter image
    if ~isempty(filterstruct) && ~isempty(fieldnames(filterstruct))
        image = tom_apply_filter(image,filterstruct);
    end
    
    mask_image=(zeros(size(image),'single'));
    mask_image(sx_ref_vol_2-sx_ref_2+1:sx_ref_vol_2+sx_ref_2, sy_ref_vol_2-sy_ref_2+1:sy_ref_vol_2+sy_ref_2)=mask_template;

    image2_fft=fft2(image.*image);
    image_fft=fft2(image);
    mask_image_fft=fft2(mask_image);

    
    T_M=(1./P).*ifftshift(real(ifft2(mask_image_fft.*image_fft))); % see equation 8,9 in Roseman 2003, Ultramicroscopy
    
    sigma2_MT=(1./P).*(ifftshift(real(ifft2(mask_image_fft.*image2_fft)))-real(T_M).^2);
    
    %preallocate output images
    angles_out = zeros(size(image,1),size(image,2),'uint16');
    template_out = zeros(size(image,1),size(image,2),'uint8');
    ccf_out = zeros(size(image,1),size(image,2),'single');
    
    %loop over angles and all templates
	lauf=1;
    for templ = 1:templatesize(3)
        for ang = angs
            tempstack_img = tom_emreadc([fftfilename '_img'],'subregion',[1 1 lauf],[size(image,1)-1 size(image,2)-1 0]);
            tempstack_real = tom_emreadc([fftfilename '_real'],'subregion',[1 1 lauf],[size(image,1)-1 size(image,2)-1 0]);
            tempstack = single(tempstack_real.Value + i .* tempstack_img.Value);
            clear('tempstack_img');
            clear('tempstack_real');
            correlation_function=ifftshift(real(ifft2(tempstack.*image_fft))); % see equation 2 in Roseman 2003, Ultramicroscopy 
            correlation_function=1./(P.*abs(sqrt(sigma2_MT))).*correlation_function; % see equation 6 in Roseman 2003, Ultramicroscopy 
            [ccf_out,indices] = max(cat(3,ccf_out,correlation_function),[],3); %get maximum of correlation function
            inds = find(indices==2);
            angles_out(inds) = uint16(ang);
            template_out(inds) = uint8(templ);
            lauf = lauf+1;
        end
    end
    
    %mask ccfmap
    ccf_out(1:sx_ref_2,:) = 0;
    ccf_out(size(ccf_out,1)-sx_ref_2:size(ccf_out,1),:) = 0;
    ccf_out(:,1:sy_ref_2) = 0;
    ccf_out(:,size(ccf_out,2)-sy_ref_2:size(ccf_out,2)) = 0;
    
    if size(maskcell,2) >= l && ~isempty(maskcell(l))
        %generate polygon mask
        mask_st.Apply=1;
        mask_st.Value.size_x=size(image,1);
        mask_st.Value.size_y= size(image,2);
        mask_st.Polygons = maskcell{l};
        mask_st.Type='roipoly';
        mask_st.Inverse=0;
        mask_st.Value.angle_phi = 0;
        mask_st.Value.angle_psi = 0;
        mask_st.Value.angle_theta = 0;
        mask_st.Value.binning = binning;
        ccfmask=tom_create_mask(mask_st);
        ccf_out = ccf_out .* imerode(ccfmask,se);
        clear ccfmask;
    end
    
    
    %extract particles
    align=tom_av2_extractparticles(ccf_out, angles_out, template_out, templatesize(1), number, filename{1}, binning, align);
    if ccfmapflag == 1
        ccfmap = ccf_out;
    else
        ccfmap = [];
    end
    disp([num2str(l) ' of ' num2str(nofiles) ' images done']);
    l=l+1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [normed_vol mea st n]=norm_inside_mask(vol,mask)

% set the mean to zero and the standard deviation to 1
%
mask=mask~=0;
n=sum(sum(mask~=0));
mea=sum(sum((vol.*mask).*(vol~=0)))./n;
st=sqrt(sum(sum((((mask==0).*mea)+(vol.*mask) -mea).^2))./n);
normed_vol=((vol-mea)./st).*mask;


function ind = segment(im)

im2 = tom_filter(im,10,'circ','real');
mask = EMSeg(im2,2);
se2 = strel('disk',7);
se = strel('disk',10);
mask2 = imdilate(imdilate(imdilate(imdilate(imerode(mask,se),se),se2),se2),se2);
%figure;tom_imagesc(mask);
ind = find(mask2 == 1);
%im(ind) = 0;
%figure;tom_imagesc(im);

