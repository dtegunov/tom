function [markerset, matchchain, comparelist] = tom_mark_findmatches(filenames, markerset, ncorrelate, comparelist, matchparam, imreadbinning, improcparam, verbose, hwaitbar)
%TOM_MARK_FINDMATCHES try to find markers in a stack of images.
% 
%   [markerset, matchchain, comparelist] = 
%                  TOM_MARK_FINDMATCHES(filenames, markerset, ncorrelate, comparelist,  
%                                       matchparam, imreadbinning, improcparam)
%   Using TOM_MARK_MATCHNCC it tries to fill up a not full
%   markerset with matched imagecoordinates. 
%   It does so, by compareing neighbouring images (from the stack).
%
%PARAMETERS
%  INPUT
%     filenames: Are the filenames of the images. It must be a cell-array
%         of strings; i.e. the filenames of the em-files. The order of the
%         images is here important to know the order of the image stack. The
%         tiltangle from the emfiles or the filename is ignored; the ordering
%         only depends on the appearence of the filenames in this
%         parameter.
%     markerset: is a matrix containg the pixelcoodinates of every marker
%         for every image from the stack. It has the dimension
%         (2 x nImages x nMarkers). It means, that the marker #m has in
%         image #i the pixelcoodinates [x y] = markerset([1 2], i, m).
%         If at least one of the two coodinates from a marker (in one
%         image) is <0 or NaN, it means, that the marker is not found in
%         this image, and is thereby to be searched by this function.
%         Even if the images are binned or resized these coordinates of the
%         markerset are always the real position to the marker in the
%         original image/emfile.
%     ncorrelate: number of neighbouring images, in which a missing
%         markerpoint is to be searched using the normalized cross correlation.
%         e.g. if the markerpoint #m in image #i is unknown, but
%         it is known in its "ncorrelate" closest neighbours, it is
%         searched there. 'Neighbours' always in regard to the ordering in the 
%         input-stack, the images.
%         The search is done beginning with the closest until a positive
%         match were found.
%         If a NCC returns successfully, the found position of the marker
%         is affilated in the markerset, and following searches tread them
%         as any other point in the markerset.
%     comparelist: Is a array of size (N x 3) where every row is a vector of
%         indeces [idx_mark, idx_imfound, idx_im_empty]. This is a list of
%         marker-image-combinations which are not to be checked. This
%         means, that the missing markerposition from marker #idx_mark in
%         the image #idx_im_empty is not searched in the image
%         #idx_imfound. 
%         The use of this field is the following: If after a run of the
%         function a false match is returned by comparing two images, 
%         the false position can be deleted and called the function again
%         (with the same matching-parameters). To prevent to rematch again 
%         the same positions, the tuppel of marker and images can be entered in
%         "comparelist".%        
%     matchparam: A structure containing the parameters passed to the
%         correlation function TOM_MARK_MATCHNCC. "matchparam" contains the
%         fields "w", "max_disparity", "min_val", "mutual" and "im_shift".
%         All parameters have the same meaning as their appropriate
%         parameters in TOM_MARK_MATCHNCC. 
%         The only exception is the field "im_shift",
%         which instead of beeing a 2-vector with [x y] displacement, is a
%         matrix of size (nImages x nImages x 2) containg all the
%         displacements between each two images. TOM_MARK_GETAVGSHIFT can
%         be used to retain these displacements.
%         Every field can be ommited, leaving the defaultvalues as
%         suggested in TOM_MARK_MATCHNCC.
%         The size from w, im_shift, etc. is given according to the real
%         image dimension. The function adjusts the values by itself, according 
%         to the used binning/resize. To do so,
%         TOM_MARK_APPLYTRANSF is called.
%     imreadbinning: Defines the binning used when reading the images with
%         TOM_EMREADC. Obviosly this is "transformatio" applied to every
%         image before doing what is defined in "improcparam".
%     improcparam: is a field of structures, containing the parameters of
%         a pre-processing of the images. After reading the image with TOM_EMREADC 
%         (binning "imreadbinning") the image is handles the the function
%         TOM_MARK_APPLYTRANSF with the structure "improcparam".
%  OUTPUT:
%     markerset: The old markerset, fill up with any further found
%         markerpoint.
%     matchchain: For every markerposition it has the index of
%         the image from were it was matched during the computation.
%         It is NaN if the marker were not found for the corresponding image
%         or 0, if the position were already know at the beginning of the
%         call.
%         Is usefull for false matches. If there is found a false markerpoint 
%         by the function, it has to be deleted/changed. As every further
%         markerposition, which derives from a match to this wrong
%         position.
%     comparelist: is the same as the inputparameter "comparelist", enhanced by
%         every done comparison during the computation, in the order as it
%         were done.
%
%   See also TOM_MARK_MATCHNCC, TOM_MARK_APPLYTRANSF,
%   TOM_MARK_GETAVGSHIFT
%
%REFERENCES
%
%   20/03/07
%
%   created by Thomas Haller
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Bring the markerset-argument in the right form.
if (isstr(markerset))
    markerset = tom_emread(markerset);
end;
if (isstruct(markerset))
    markerset = markerset.Value;
end;
msize = int32(size(markerset));
if (length(msize)==2)
    msize(3) = 1;
end;
if (length(filenames) ~= msize(2))
    error(['markerset and filenames have not the same number of images.']);
end;
if (msize(1) ~= 2)
    markerset = markerset([2 3], :, :);
    msize(1) = 2;
end;



% check input parameters and/or set defaults...
if (~exist('verbose', 'var') || ~islogical(verbose) || numel(verbose)~=1)
    verbose = true;
end;
if (~exist('hwaitbar','var') || ~ishandle(hwaitbar) || numel(hwaitbar)~=1)
    hwaitbar = [];
end;
if (isempty(imreadbinning) || isnan(imreadbinning) || uint8(imreadbinning)~=imreadbinning)
    imreadbinning = 0;
end;
if (~(ncorrelate > 0))
    ncorrelate = 1;
end;
% Prepareing matchparameters and setting default values for ommited fields.
if (~isstruct(matchparam) || isempty(matchparam))
    matchparam = struct();
end;
if (~isfield(matchparam, 'w'))
    matchparam.w = 'default';
end;
if (~isfield(matchparam, 'max_disparity'))
    matchparam.max_disparity = 'default';
end;
if (~isfield(matchparam, 'im_shift') || (isstr(matchparam.im_shift) && strcmp(matchparam.im_shift, 'default')))
    matchparam.im_shift = tom_mark_getAvgShift(msize(2));
end;
if (~isfield(matchparam, 'min_val'))
    matchparam.min_val = 'default';
end;
if (~isfield(matchparam, 'mutual'))
    matchparam.mutual = 'default';
end;


if (isempty(comparelist) || size(comparelist, 2)~=3)
    comparelist = zeros(2*ncorrelate*msize(2)*msize(3), 3, 'int32');
    ncomparelist = 0;
else
    ncomparelist = size(comparelist, 1);
    comparelist = cat(1, comparelist, zeros(2*ncorrelate*msize(2)*msize(3) - ncomparelist, 3, 'int32'));
end;


ncorrelate = int32(ncorrelate);

% Read the headers of images. and prepare the image structure.
imagestruct = tom_mark_loadEMImages('init', filenames);
imagestruct.emread_binning = imreadbinning;
imagestruct.transform = improcparam;
%imagestruct.freesize = 2500000;
imagestruct = tom_mark_loadEMImages('header', imagestruct, 1:msize(2));


idx_existing_markers = all(markerset>0, 1);

nfound_markerpoints = uint32(sum(sum(idx_existing_markers)));

markerset(:, ~idx_existing_markers) = nan;

% Transform the markerset and the matchparams according the binning and
% image-processing (like resizeing).
markerset = tom_mark_applyTransf('markerset_pre', markerset, improcparam, imreadbinning);
matchparam = tom_mark_applyTransf('matchparam_pre',matchparam, improcparam, imreadbinning);


% Find for every missing markerposition, the distance (in the imagestack)
% to the closest known markerposition and insert it into the matrix
% marker_comparedistance.
marker_comparedistance = zeros(1, msize(2), msize(3), 'int32');


for (ii=1:msize(2))
    if (~imagestruct.error(ii))
        for (im=1:msize(3))
            if (~idx_existing_markers(1, ii, im))
                distance = ncorrelate+1;
                for (k=1:ncorrelate)
                    kidx_low  = ii-k;
                    kidx_high = ii+k;
                    if ((kidx_low  >  0        && idx_existing_markers(1, kidx_low , im)) || ...
                        (kidx_high <= msize(2) && idx_existing_markers(1, kidx_high, im)))
                        distance = k;
                        break;
                    end;
                end;
                marker_comparedistance(1, ii, im) = distance;
            end;
        end;
    else
        marker_comparedistance(1, ii, 1:msize(3)) = -1;
    end;
end;
%marker_comparedistance(1, ~idx_existing_markers) = 1;

%markerset(:, ~idx_existing_markers) = -1;%%%
 



matchchain = nan(1, msize(2), msize(3));
matchchain(1, idx_existing_markers) = 0;



neighbourhood = [-ncorrelate:-1, 1:ncorrelate];

ic = int32(1);

sortdirection_str = 'ascend';
sortdirection_num = [ 1 2 ];

nlastdisp = 0;

max_integerformatstring_length = fix(log10(double(msize(2)*msize(3))))+2;
max_integerformatstring = ['%' num2str(max_integerformatstring_length) 'd'];

distcnt_text{1} = ['           ' repmat(' ',1,max_integerformatstring_length) 'faild:                      '];
distcnt_text{2} = ['           ' repmat(' ',1,max_integerformatstring_length) 'successfull:                '];
distcnt_text{3} = ['           ' repmat(' ',1,max_integerformatstring_length) 'skip because already found: '];
distcnt_text{4} = ['           ' repmat(' ',1,max_integerformatstring_length) 'skip because closer found:  '];

breakloop = false;

try
    oldtoc = toc;
catch
    tic;
    oldtoc = toc;
end;

while (ic <= ncorrelate)
    
    if (verbose)
        s = '';
        s = ['%%% ' num2str(toc()) ' seconds %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%' char(10)];
        disp([repmat(char(8),1,nlastdisp) s sprintf('%7.3f', double(100*nfound_markerpoints)/double(msize(2)*msize(3))), '%: ...']);
        nlastdisp = 0;
        v_distancecnt_per_run = zeros(5, ncorrelate, 'int32');
    end;
    if (ishandle(hwaitbar))
        waitbar(double(nfound_markerpoints)/double(msize(2)*msize(3)), hwaitbar);
        if (get(hwaitbar, 'UserData'))
            breakloop = true;
            break;
        end;
    end;

    
    % Generate a list of all comparisions which have to be made.
    correlate_list = zeros(0, 3, 'int32');
    ncorrelate_list = int32(0);

    marker_comparedistance_new = marker_comparedistance;
    
    
    for (ii=1:msize(2))
        if (~imagestruct.error(ii))
            for (im=1:msize(3))
                cd = marker_comparedistance(1,ii,im);
                if (cd <= 0)
                elseif (cd == ic)
                    j = ii-ic;
                    if (j >= 1 && marker_comparedistance(1,j,im)==0)
                        ncorrelate_list = ncorrelate_list+1;
                        correlate_list(ncorrelate_list, [1 2 3])  = [im j ii];
                    end;
                    j = ii+ic;
                    if (j<=msize(2) && marker_comparedistance(1,j,im)==0)
                        ncorrelate_list = ncorrelate_list+1;
                        correlate_list(ncorrelate_list, [1 2 3])  = [im j ii];
                    end;
                    marker_comparedistance_new(1,ii,im) = cd+1;
                elseif (cd <= ncorrelate)
                    found = false;
                    for (k=1:cd-1)
                        j = ii - k;
                        if (j>0 && marker_comparedistance(1, j, im)>0 && marker_comparedistance(1, j, im)<cd)
                            found = true;
                            break;
                        end;
                        j = ii + k;
                        if (j<=msize(2) && marker_comparedistance(1, j, im)>0 && marker_comparedistance(1, j, im)<cd)
                            found = true;
                            break;
                        end;
                    end;
                    if (~found)
                        j = ii-cd;
                        if (j >= 1 && marker_comparedistance(1,j,im)==0)
                            ncorrelate_list = ncorrelate_list+1;
                            correlate_list(ncorrelate_list, [1 2 3])  = [im j ii];
                        end;
                        j = ii+cd;
                        if (j<=msize(2) && marker_comparedistance(1,j,im)==0)
                            ncorrelate_list = ncorrelate_list+1;
                            correlate_list(ncorrelate_list, [1 2 3])  = [im j ii];
                        end;
                        marker_comparedistance_new(1,ii,im) = cd+1;
                    end;
                else
                end;
            end;
        end;
    end;

    correlate_list = setdiff(correlate_list, comparelist(1:ncomparelist, :), 'rows');
    
    [image_pairs, idx_image_pairs] = sortrows(sort(correlate_list(:,2:3), 2, sortdirection_str), sortdirection_num);
    correlate_list = correlate_list(idx_image_pairs, :);
    
 
    if (strcmp(sortdirection_str, 'ascend'))
        sortdirection_str = 'descend';
        sortdirection_num = -[ 1 2 ];
    else
        sortdirection_str = 'ascend';
        sortdirection_num = [ 1 2 ];
    end;
    
    
    ncorrelate_list = int32(size(correlate_list, 1));
    

    nextic = ic+1;
    

    
    if (verbose)
        idx = abs(correlate_list(:,2) - correlate_list(:,3));
        for (i=1:size(idx,1))
            v_distancecnt_per_run(1, idx(i,1)) = v_distancecnt_per_run(1, idx(i,1))+1;
        end;
        s = '';
        for (i=1:ncorrelate)
            s = [s ' ' sprintf(max_integerformatstring, v_distancecnt_per_run(1, i))];
        end;
        disp([char([8 8 8 8]) sprintf(max_integerformatstring, ncorrelate_list) ' comparisons          dist = [' s ' ]']);
    end;


    for (i=1:ncorrelate_list)
        if (toc() - oldtoc > 1)
            oldtoc = toc();
            if (verbose)
                s = num2str(100*double(i)/double(ncorrelate_list),  '%7.3f');
                disp([repmat(char(8),1,nlastdisp) s '%']);
                nlastdisp = length(s)+2;
            end;
            if (ishandle(hwaitbar))
                waitbar(double(nfound_markerpoints)/double(msize(2)*msize(3)), hwaitbar);
                breakloop = get(hwaitbar, 'UserData');
                if (breakloop)
                    break;
                end;
            end;
        end;
        c = correlate_list(i,:);

        if (marker_comparedistance_new(1,c(3),c(1)) >= marker_comparedistance(1,c(3),c(1)))

            if (isempty(imagestruct.Value{c(2)}))
                imagestruct = tom_mark_loadEMImages(imagestruct, c(2));
            end;
            if (isempty(imagestruct.Value{c(3)}))
                imagestruct = tom_mark_loadEMImages(imagestruct, c(3));
            end;
            

            [mpt, mval] = tom_mark_matchncc(imagestruct.Value{c(2)}, imagestruct.Value{c(3)}, markerset(:,c(2),c(1))', ...
                                            matchparam.w, matchparam.max_disparity, squeeze(matchparam.im_shift(c(2), c(3), :))', matchparam.min_val, matchparam.mutual);
            ncomparelist = ncomparelist + 1;
            comparelist(ncomparelist, 1:3) = c;

            
            if (isempty(mpt))
                if (verbose)
                    %disp([num2str(i) ': Search marker #' num2str(c(1)) ' from image #' num2str(c(2)) ' in #' num2str(c(3)) '... nothing found ***********************************']);nlastdisp = 0;
                    v_distancecnt_per_run(2, abs(c(3)-c(2))) = v_distancecnt_per_run(2, abs(c(3)-c(2))) + 1;
                end;
            else
                if (verbose)
                    %disp([num2str(i) ': Search marker #' num2str(c(1)) ' from image #' num2str(c(2)) ' in #' num2str(c(3)) '... marker found ####################################']);nlastdisp = 0;
                    v_distancecnt_per_run(3, abs(c(3)-c(2))) = v_distancecnt_per_run(3, abs(c(3)-c(2))) + 1;
                end;

                markerset([1 2],c(3),c(1)) = mpt';                       
                marker_comparedistance_new(1,c(3),c(1)) = 0;                       
                matchchain(1,c(3),c(1)) = c(2);
                nfound_markerpoints = nfound_markerpoints+1;
                
                % Adjust the marker_comparedistance for the neighouring
                % markerpoints...
                % First down the stack...
                cc = c(3) - 1;
                ii = 1;
                while (ii <= ncorrelate)
                    if (cc <= 0) 
                        break;
                    end;
                    if (marker_comparedistance_new(1, cc, c(1)) > ii)
                        marker_comparedistance_new(1, cc, c(1)) = ii;
                        if (nextic > ii)
                            nextic = ii;
                        end;
                        ii = 1;
                    else
                        ii = ii + 1;
                    end;                            
                    cc = cc - 1;
                end;
                % then up...
                cc = c(3) + 1;
                ii = 1;
                while (ii <= ncorrelate)
                    if (cc > msize(2)) 
                        break;
                    end;
                    if (marker_comparedistance_new(1, cc, c(1)) > ii)
                        marker_comparedistance_new(1, cc, c(1)) = ii;
                        if (nextic > ii)
                            nextic = ii;
                        end;
                        ii = 1;
                    else
                        ii = ii + 1;
                    end;                            
                    cc = cc + 1;
                end;
            end;
        else
            if (verbose)
                if (marker_comparedistance_new(1,c(3),c(1)) == 0)
                    %disp([num2str(i) ': Search marker #' num2str(c(1)) ' from image #' num2str(c(2)) ' in #' num2str(c(3)) ' skipped this marker were already found.']);nlastdisp = 0;
                    v_distancecnt_per_run(4, abs(c(3)-c(2))) = v_distancecnt_per_run(4, abs(c(3)-c(2))) + 1;
                else
                    %disp([num2str(i) ': Search marker #' num2str(c(1)) ' from image #' num2str(c(2)) ' in #' num2str(c(3)) ' skipped because there was found a closer marker to check first.']);nlastdisp = 0;
                    v_distancecnt_per_run(5, abs(c(3)-c(2))) = v_distancecnt_per_run(5, abs(c(3)-c(2))) + 1;
                end;
            end;
        end;                
    end;
    
    if (breakloop)
        break;
    end;
    
    if (verbose)
        s = '';
        for (i=1:4)
            s = [s distcnt_text{i} '['];
            for (j=1:ncorrelate)
                s = [s ' ' sprintf(max_integerformatstring, v_distancecnt_per_run(i+1, j))];
            end;
            s = [s ' ]' char(10)];
        end;
        disp([repmat(char(8),1,nlastdisp) s(1:end-1)]);
        nlastdisp = 0;

        if (~breakloop && any(v_distancecnt_per_run(1,:) ~= sum(v_distancecnt_per_run(2:5,:))))
            warning('UNEXPECTED');
        end;
    end;
    

    marker_comparedistance = marker_comparedistance_new;
     
    ic = nextic;

%    if (~isempty(marker_comparedistance_new(marker_comparedistance_new(1,:,:) > 0)) && ic > min(marker_comparedistance_new(marker_comparedistance_new(1,:,:) > 0)))
%        error('UNEXPECTED POSITION IN CODE');
%    end;
end;
if (verbose)
    disp([repmat(char(8),1,nlastdisp)]);
end;



if (any(~isnan(markerset([1 2], marker_comparedistance_new~=0))))
    error('UNEXPECTED POSITION IN CODE');
    %markerset([1 2], marker_comparedistance_new~=0) = nan;
end;

comparelist = comparelist(1:ncomparelist, 1:3);
markerset = tom_mark_applyTransf('markerset_post', markerset, improcparam, imreadbinning);





