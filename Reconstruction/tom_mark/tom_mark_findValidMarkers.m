function tom_mark_findValidMarkers(markerset)


msize = size(markerset);
if (length(msize) < 3)
    msize(3) = 1;
end;



marker_exists = all(markerset > 0, 1);


for (i=1:10)

    
    idx = extract_subset(marker_exists, 10, 5);

%    mf = 
% TODO   
%    [Matrixmark, psi, sigma, x, y, z]  = tom_alignment3d(Matrixmark, irefmark, imintilt, r, imdim)
    
end;



return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Returns a index vector to the 3 dimenstion of the markerset.
% First The index is chosen randomly, so that at least 
% min(msize(3), n) elements are returned.
% It its possible then the result is augmented until for
% every image at least minnimage markerpoints are in the 
% result.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function idx = extract_subset(marker_exists, n, minnimage)

msize = size(marker_exists);
if (msize(3) <= n)
    idx = 1:msize(3);
    return;
end;


idx = randperm(msize(3));
idx = idx(1:n);
idxnotmask = true(1,msize(3));
idxnotmask(idx) = false;


imageidx = randperm(msize(2));
if (minnimage > 0)
    for (i=imageidx)
        nmissing = minnimage - sum(marker_exists(1,i,idx));
        if (nmissing > 0)
            idxcandidates = find(and(squeeze(marker_exists(1,i,:))', idxnotmask));
            if (nmissing >= length(idxcandidates))
            else
                foundidxrand = randperm(length(idxcandidates));
                idxcandidates = idxcandidates(foundidxrand(1:nmissing));
            end;            
            idx = [idx idxcandidates];
            idxnotmask(idxcandidates) = false;
        end;
    end;
end;



