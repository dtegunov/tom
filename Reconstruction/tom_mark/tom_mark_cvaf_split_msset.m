function markerset = tom_mark_cvaf_split_msset(markerset, threshold, nTrails, hwaitbar)
% Splits a markerset.

markerset(3,:,:) = 1;

msize = size(markerset);
if (length(msize)==2)
    msize(3) = 1;
end;

if (~exist('hwaitbar','var') || ~ishandle(hwaitbar) || numel(hwaitbar)~=1)
    hwaitbar = [];
end;


maxlength = uint16(repmat(msize(2):-1:1, [1, 1, msize(3)]));
for (i=1:msize(3))
    base = 0;
    for (j=1:msize(2))
        if (any(isnan(markerset(:,j,i)) | isinf(markerset(:,j,i))))
            maxlength(1, (base+1):j, i) = ((j-1):-1:base)-base;
            base = j;
        end;
    end;
end;

tic;
for (setn = 2:msize(2))
    
    if (ishandle(hwaitbar))
        waitbar(setn / msize(2), hwaitbar);
        if (get(hwaitbar, 'UserData'))
            markerset = [];
            return;
        end;
    end;
    
    
    for (iim = 0:(msize(2)-setn))
        iims = iim+(1:setn);
        
        relevant_idxm = squeeze(maxlength(1, iim+1, :) >= setn)';
        relevant_idx = find(relevant_idxm);
        if (~isempty(relevant_idx))
            ms = markerset(:, iims, relevant_idxm);

            [idxinliers] = tom_mark_cvaf_inliersRANSAC(ms, threshold, nTrails, true);

            idxinliers = relevant_idx(idxinliers);
            idxoutliers = setdiff(relevant_idx, idxinliers);

            for (i=idxoutliers)
                startidx = iim+1;
                while (startidx >= 1 && maxlength(1, startidx, i) > setn+iim-startidx)                
                    startidx = startidx - 1;
                end;
                maxlength(1, (startidx:iim)+1, i) = (iim:-1:startidx) - startidx + setn - 1;
            end;
        end;
    end;
end;

maxl_peaks = maxlength(1,2:msize(2), :);
maxl_peaks(maxlength(1,1:(msize(2)-1),:) > maxl_peaks) = 0;
maxl_peaks = cat(2, maxlength(1,1,:), maxl_peaks);
maxl_peaks(maxl_peaks == 1) = 0;

newms = nan(2, msize(2), length(find(maxl_peaks > 0)));
msize3 = 0;
for (i=1:msize(3))
    maxl_peak = uint16(find(maxl_peaks(1,:,i)));
    [j, k] = sort(maxl_peaks(1,maxl_peak, i), 'descend');
    maxl_peak = maxl_peak(k);
    idxs = false(1, msize(2), length(maxl_peak));
    for (j=1:length(maxl_peak))
        idxs(1, maxl_peak(j)-1+(1:maxl_peaks(1,maxl_peak(j),i)), j) = true;
    end;
    
    % If there are intersections of chains take the longest ones...
    % (maybe it would be better to take the "best" soltution where in total
    % the longes chains are, witch taking no point more then often. but
    % difficult to implment....
    % Or just take all the chains? and count (i.e. weight) some markers
    % more????
    for (j=2:length(maxl_peak))
        idxs(1,find(sum(idxs(1,:,1:j), 3) > 1), j) = 0;
    end;
    
   idxs = idxs(1, :, sum(idxs, 2) >= 2);
   
    for (j=1:size(idxs,3))
        msize3 = msize3 + 1;
        newms(1:2, idxs(1,:,j), msize3) = markerset(1:2,idxs(1,:,j), i);
    end;    
end;

markerset = newms(1:2, :, 1:msize3);


