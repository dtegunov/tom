function [markerset, match_chain, nocompare] = tom_mark_delMarkerPos(dellist, markerset, match_chain, nocompare)
%TOM_MARK_DELMARKERPOS deletes markerpositions from a markerset.
% 
%   [markerset, match_chain, nocompare] = 
%        tom_mark_delMarkerPos(dellist, markerset, match_chain, nocompare)
%
%   Deletes from a markerset the position from the marker in the image and
%   also the positions which derive from a match to the deleted one (i.e. a
%   cascated deleting).
%
%PARAMETERS
%  INPUT
%     dellist: Is a matrix of indexpairs of size N x 2. It specifies N
%         positions which are to delete, where the first row is the index of
%         the marker in the markerset, and the second is the index of the
%         image.
%     markerset: a matrix of size 2 x nImages x nMarker where the
%         pixelcoodinates of each marker in each image is saved.
%         Deleting a marker in an image, means to set the coordinates to
%         NaN.
%     match_chain: A matix of size 1 x nImages x nMarker as returned from
%         TOM_MARK_FINDMATCHES
%     nocompare: A matrix of size M x 3.
%  OUTPUT:
%     markerset: The inputargument, where the positions named in @dellist@
%         are set to NaN. Further there are all posistions set to NaN,
%         which derive from the deleted positions, according to
%         "match_chain"
%     match_chain: The input-argument where the deleted positions are set
%         to NaN.
%     nocompare: The input argument reduced of the list of matches made
%         agains deleted markerpositions.
%
%   See also TOM_MARK_FINDMATCHES
%
%REFERENCES
%
%   23/03/07
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


msize = int32(size(markerset));
if (length(msize) == 2)
    msize(3) = 1;
end;

if (msize(2) ~= size(match_chain, 2) || msize(3) ~= size(match_chain, 3))
    error('Dimension of the markerset and the match_chain does not match.');
end;

if (isempty(dellist))
    return;
end;

if (size(dellist, 2) ~= 2)
    error('The list of points to delete is not Nx2. Every row has the meaning [MarkerIndex x ImageIndex]');
end;

ndellist = size(dellist, 1);

idx_dellist = all([all(dellist > 0, 2), dellist(:,1)<=msize(3), dellist(:,2)<=msize(2)], 2);
nocompare_idx = true(size(nocompare, 1),1);
for i=1:ndellist
    if (idx_dellist(i))
        im = dellist(i, 1);
        ii = dellist(i, 2);
        
        todelete = [ ii ];
        ntodelete = 1;

        while (~isempty(todelete))
            k = todelete(1);
            todelete(1) = [];

            markerset([1 2], k, im) = nan;
            match_chain(1, k, im) = nan;
            if (k==ii)
                nocompare_idx = nocompare_idx & (nocompare(:,1)~=im | (nocompare(:,2)~=k                    ));
            else
                nocompare_idx = nocompare_idx & (nocompare(:,1)~=im | (nocompare(:,2)~=k & nocompare(:,3)~=k));
            end;

            todelete = [todelete find(match_chain(1, 1:msize(2), im) == k)];
            
        end;
    end;
end;

nocompare = nocompare(nocompare_idx, [1 2 3]);
