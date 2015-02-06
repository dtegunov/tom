function [markerset_new] = tom_mark_refineMarkerSet(filenames, markerset, ncorrelate, matchparam, imreadbinning, improcparam)

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



% check input parameters and set defaults...
if (~isnumeric(imreadbinning) || numel(imreadbinning)~=1 || uint8(imreadbinning)~=imreadbinning)
    imreadbinning = 0;
end;
if (~isnumeric(ncorrelate) || uint8(ncorrelate)~=ncorrelate)
    ncorrelate = 1;
end;
if (numel(ncorrelate) == 1)
    neighbourhood = [-ncorrelate:-1, 1:ncorrelate];
else
    neighbourhood = ncorrelate;
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
if (~isfield(matchparam, 'min_val'))
    matchparam.min_val = 'default';
end;
if (~isfield(matchparam, 'mutual'))
    matchparam.mutual = 'default';
end;
if (~isfield(matchparam, 'im_shift'))
    matchparam.im_shift = 'default';
end;

ncorrelate = int32(ncorrelate);

% Read the headers of images.
imagestruct = tom_mark_loadEMImages('init', filenames);
imagestruct.emread_binning = imreadbinning;
imagestruct.transform = improcparam;
%imagestruct.freesize = 2500000;
imagestruct = tom_mark_loadEMImages('header', imagestruct, 1:msize(2));


idx_existing_markers = all(markerset>0, 1);

nfound_markerpoints = uint32(sum(sum(idx_existing_markers)));
nfound_markerpoints_all = uint32(msize(2)*msize(3));


markerset_new = zeros(2, msize(2), msize(3));
markerset_new_cnt = zeros(1, msize(2), msize(3));


% Transform the markerpositions accoding to the binning and
% image-processing (like resizeing).
markerset = tom_mark_applyTransf('markerset_pre', markerset, improcparam, imreadbinning);
matchparam = tom_mark_applyTransf('matchparam_pre', matchparam, improcparam, imreadbinning);


nlastdisp = 0;

for (ii=1:msize(2))
    
    s = sprintf('%7.3f%%', 100*double(ii)/double(msize(2)));
    disp([repmat(char(8),1,nlastdisp) s]);
    nlastdisp = length(s)+1;

    if (imagestruct.error(ii))
        continue;
    end; 
    if (isempty(imagestruct.Value{ii}))
        imagestruct = tom_mark_loadEMImages(imagestruct, ii);
    end;
    
    for (j=neighbourhood)
        ji = ii+j;
        if (ji<1 || ji>msize(2) || imagestruct.error(ji))
            continue;
        end;
        if (isempty(imagestruct.Value{ji}))
            imagestruct = tom_mark_loadEMImages(imagestruct, ji);
        end;
        for (im=1:msize(3))
            if (idx_existing_markers(1,ii,im) && idx_existing_markers(1,ji,im))
                [mpt, mval] = tom_mark_matchncc(imagestruct.Value{ji}, imagestruct.Value{ii}, markerset(:,ji,im)', ...
                                            matchparam.w, matchparam.max_disparity, (markerset(:,ii,im) - markerset(:,ji,im))', matchparam.min_val, matchparam.mutual);
        
                if (~isempty(mpt))
                    markerset_new([1 2], ii, im) = markerset_new([1 2], ii, im) + mpt';
                    markerset_new_cnt(1, ii, im) = markerset_new_cnt(1, ii, im) + 1;
                end;
            end;        
        end;
    end;
end;

markerset_new(:,markerset_new_cnt <= 0) = nan;
markerset_new(1,:,:) = markerset_new(1,:,:) ./ markerset_new_cnt(1,:,:);
markerset_new(2,:,:) = markerset_new(2,:,:) ./ markerset_new_cnt(1,:,:);

markerset_new = tom_mark_applyTransf('markerset_post', markerset_new, improcparam, imreadbinning);





