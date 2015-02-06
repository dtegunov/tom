function [new_markerset, match_chain_new] = tom_mark_splitmarkerset(markerset, nsplit, minnfound, match_chain, overlap)

    
msize = size(markerset);
if (length(msize) == 2)
    msize(3) = 1;
end;


if (~exist('nsplit', 'var') || prod(size(nsplit)) ~= 1 || nsplit~=round(nsplit) || ~(nsplit < msize(2)/2 && nsplit >= 1))
    nsplit = floor(msize(2) / 5);
end;


lsplit = double(msize(2)-1) / double(nsplit);


idx_low = 1;
idx = zeros(nsplit, 2);
for (i=1:nsplit)
    idx_high = floor((i)*lsplit) + 1;
    
    idx(i, 1) = idx_low;
    idx(i, 2) = idx_high;

    idx_low = idx_high + 1;
end;

if (exist('overlap', 'var') && overlap)
    idx(2:nsplit  ,1) = idx(2:nsplit  ,1) - 1;
    idx(1:nsplit-1,2) = idx(1:nsplit-1,2) + 1;
end;


new_markerset = nan(2, msize(2), msize(3)*nsplit);

insplit = 1;
for (i=1:msize(3))
    marker = markerset([1 2], :, i);
    for (j=1:nsplit)
        idxline = idx(j,1) : idx(j,2);
        new_markerset([1 2], idxline, insplit) = marker([1 2], idxline, 1);
        insplit = insplit + 1;
    end;
end;    


if (exist('match_chain_new', 'var') && ~isempty(match_chain_new))
    match_chain_new = nan(1, msize(2), msize(3)*nsplit);
    insplit = 1;
    for (i=1:msize(3))
        marker = markerset([1 2], :, i);
        match_chain_row = match_chain(1, :, i);
        for (j=1:nsplit)
            idxline = idx(j,1) : idx(j,2);
            match_chain_idx = any([match_chain_row(1,idxline,1)==0; ismember(match_chain_row(1,idxline,1), idxline)],1);
            match_chain_idx0 = (~match_chain_idx) & all(marker([1 2], idxline, 1) > 0, 1);
            match_chain_new(1,idxline(match_chain_idx0),insplit) = 0;
            match_chain_new(1,idxline(match_chain_idx),insplit) = match_chain_row(1,idxline(match_chain_idx),1);
            insplit = insplit + 1;
        end;
    end;    
else
    match_chain_new = [];
end;


if (~exist('minnfound', 'var') || prod(size(minnfound)) == 1 && minnfound >= 1)
    idx = find(sum(all(new_markerset > 0, 1)) >= minnfound);
    new_markerset = new_markerset([1 2], 1:msize(2), idx);
    
    if (~isempty(match_chain_new))
        match_chain_new = match_chain_new(1, 1:msize(2), idx);
    end;
    
end;







