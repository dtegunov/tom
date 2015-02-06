function [P, X, aligncfg] = tom_mark_cvaf_affineTraingulation(markerset, aligncfg, abort_factor, maxloop, P0, hwaitbar)

if (~exist('abort_factor', 'var') || numel(abort_factor) ~= 1)
    abort_factor = 5;
end;
if (~exist('hwaitbar','var') || numel(hwaitbar)~=1 || ~ishandle(hwaitbar) )
    hwaitbar = [];
end;
if (~exist('maxloop') || numel(maxloop)~=1)
    maxloop = 100;
end;

ni = size(markerset, 2);
msize3 = size(markerset, 3);


% Initialise random cameras and triangulate the points for the first time.
if (~exist('P0','var') || isempty(P0) || size(P0,1)~=3 || size(P0,2)~=4 || size(P0,3) ~= ni)
    P = rand(2,4,ni); 
else
    P = P0;
end;
P(3,1:3,:) = 0; P(3,4,:) =1;
X = tom_mark_cvaf_triX(P, markerset); 



mindist = Inf;
count = 0;

x = nan(3, ni, msize3);

count_loop = 0;

while (count < 20 && count_loop < maxloop)
    count_loop = count_loop + 1;
    
    P = tom_mark_cvaf_estimateP(X, markerset);
    X = tom_mark_cvaf_triX(P, markerset);
    

    % Reproject every 3D point with the computed cameras...
    for (i=1:ni)
        x(1:3, i, :) = P(:,:,i) * X;
    end;
    

    % Calculate geometric distance from original points 
    % (i.e. reprojection error)
    dist = sqrt(sum((x(1:2,:,:) - markerset) .^ 2, 1));
    ndist = sum(sum(dist(isfinite(dist)))) /  sum(sum(isfinite(dist)));


    s = [num2str(count_loop) ': meandist=' sprintf('%10.4f', ndist)];
    if (ndist < mindist)
        Pmin = P;
        Xmin = X;
        
        improovement = 100*(1-(ndist / mindist));
        count = max(0, count - abort_factor*improovement);
        mindist = ndist;
        s = [s ' improvement=' sprintf('%7.4f', improovement) '%'];
    end;
    disp(s);

    count = count + 1;
    
    if (ishandle(hwaitbar))
        waitbar(count / maxloop, hwaitbar);
        if (get(hwaitbar, 'UserData'))
            break;
        end;
    end;
end;

X = Xmin;
P = Pmin;
aligncfg =[];








