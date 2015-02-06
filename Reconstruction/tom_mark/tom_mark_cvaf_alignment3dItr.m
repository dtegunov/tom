function [P, X, aligncfg] = tom_mark_cvaf_alignment3dItr(markerset, aligncfg, abort_factor, maxloop, initshift, hwaitbar)

msize = zeros(1, 3);
[msize(1), msize(2), msize(3)] = size(markerset);

tic();

if (ischar(aligncfg) && strcmp(aligncfg, 'meanshift'))
    P = getMeanShift(markerset);
    return;
end


if (~exist('abort_factor', 'var') || numel(abort_factor) ~= 1)
    abort_factor = 5;
end;

if (~exist('hwaitbar','var') || numel(hwaitbar)~=1 || ~ishandle(hwaitbar) )
    hwaitbar = [];
end;

if (~exist('maxloop','var') || numel(maxloop)~=1)
    maxloop = 50;
end;

if (~exist('initshift','var') || size(initshift,1)~=2 || size(initshift,2)~=msize(2) || size(initshift,3)~=1)
    initshift = getMeanShift(markerset);
    initshift = initshift - repmat(initshift(:, find(abs(aligncfg.tiltangles) == min(abs(aligncfg.tiltangles)), 1), 1), [1, msize(2), 1]) + aligncfg.imsize/2;
end;

aligncfg.psi = 0;
aligncfg.tx = initshift(1,:);
aligncfg.ty = initshift(2,:);





mindist = inf;
count = 0;
count_loop = 0;
 
s = repmat(' ', [1, 256]);

irefX = msize(3)+1;
x_proj = nan(2, msize(2), irefX);
x_proj(1:2,:,irefX) = initshift;

rigAlignmin = aligncfg;
Pmin = nan(3,4,msize(2));


while (count < 12000 && count_loop < maxloop)
    count_loop = count_loop + 1;
    
    markerset(1:2, :, irefX) = x_proj(1:2,:,irefX);
    

    % The default 3D coordinate of the reference marker.

    % do the alingment...
%    [Matrixmark, align3d_psi, align3d_sigma, align3d_x, align3d_y, align3d_z]  = tom_alignment3d(Matrixmark, irefmark, imintilt, irefX, conf.imsize);
    [aligncfg.psi, aligncfg.tx, aligncfg.ty, aligncfg.sigma, X]  = ...
        tom_mark_cvaf_alignment3d(markerset, aligncfg.tiltangles, irefX, aligncfg.imintilt, aligncfg.imsize);
    
    X(1:3, irefX) = mean(X(:,all(isfinite(X(:,1:(irefX-1))),1)), 2);
    
    [P, x_proj] = tom_mark_cvaf_alignment3d_reproj(aligncfg.tiltangles, aligncfg.psi, aligncfg.tx, aligncfg.ty, aligncfg.imsize, X);
    
    figure(1);
    hold on;
    plot(aligncfg.tx);
    figure(2);
    hold on;
    plot(aligncfg.ty);
    
    
    % Calculate geometric distance from original points 
    % (i.e. reprojection error)
    distance = sqrt(sum((x_proj(:,:,1:end-1) - markerset(:,:,1:end-1)) .^ 2, 1));
    ndist = mean(distance(isfinite(distance)));
    

    %improovement = 100*(1-(ndist / mindist));
    improovement = mindist - ndist;
    if (ndist < mindist)
        rigAlignmin = aligncfg;
        Pmin = P;
        
        count = max(0, count - abort_factor*improovement);
        mindist = ndist;
    else
%        x = xmin;        
%        P = Pmin;
%        aligncfg = rigAlignmin;
    end;
    s = [sprintf('%4d', count_loop) ': ' sprintf('%6.2f', toc) ' sec; meandist=' sprintf('%20.16f', ndist) '; psi=' sprintf('%12.7f', aligncfg.psi*180/pi) '; c' num2str(count) '; improvement=' sprintf('%20.17f', improovement) ''];
    disp(s);

    count = count + 1;

    
    
    if (ishandle(hwaitbar))
        waitbar(count / maxloop, hwaitbar);
        if (get(hwaitbar, 'UserData'))
            break;
        end;
    end;
    
end;





aligncfg = rigAlignmin;
P = Pmin;
[Pnn, xdummy, X, xnn] = tom_mark_cvaf_alignment3d_reproj(aligncfg.tiltangles, aligncfg.psi, aligncfg.tx, aligncfg.ty, aligncfg.imsize, [], markerset(:,:,1:msize(3)), false);
X = X(:,1:end);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function shifts = getMeanShift(markerset)

msize = nan(1,3);
[msize(1), msize(2), msize(3)] = size(markerset);

markerset_shift = cat(2, nan(2,1,msize(3)), markerset) - cat(2, markerset, nan(2,1,msize(3)));
shifts = zeros(2,msize(2),1);

markerset_shifti_mean = shifts(1:2,1,1);
for (i=2:msize(2))
    markerset_shifti = markerset_shift(:,i,:);
    markerset_shifti = markerset_shifti(:,:,all(isfinite(markerset_shifti),1));
    if (~isempty(markerset_shifti))
        markerset_shifti_mean = markerset_shifti_mean - mean(markerset_shifti, 3);
    end;
    shifts(1:2,i,1) = markerset_shifti_mean;    
end;




