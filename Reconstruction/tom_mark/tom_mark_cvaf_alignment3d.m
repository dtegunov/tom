function [psi, tx, ty, sigma, X]  = tom_mark_cvaf_alignment3d(markerset, tiltangles, irefmark, irefmarkX, imdim)
%TOM_MARK_CVAF_ALIGNMENT3D: computes the rigid body alignment from a
% markerset, in a very similar way to TOM_ALIGNMENT3D 
%
% [psi, tx, ty, sigma, X]  = 
%    tom_mark_cvaf_alignment3d(markerset, tiltangles, irefmark, irefmarkX, imdim);
%
% TOM_MARK_CVAF_ALIGNMENT3D calculates the tilt-axis, the shifts and fit-errors of an
% EM markerfile
%
%PARAMETERS
%
%  INPUT
%   markerset:          A 2xNxM matrix of marker-coordinates. N is the
%                       number of projections (tilts) and M the number of
%                       marks. Unknown positions are set to Nan! (in
%                       contrast to TOM_ALIGNMENT3D where the unknown
%                       positions are set to -1.
%   tiltangles          A vector of all N tiltangles in Degrees.
%   irefmark            Index of the reference marker in the markerset (1...M)
%   irefmarkX           3D point of the reference marker i.e. a 3-vector.
%                       or:
%                       - Set irefmarkX to an integer number 1...N to use the
%                         default 3D coordinate
%                         [markerset(1:2,irefmarkX,irefmark); imdim/2 + 1].
%                       - Set irefmarkX to [] to choose as reference
%                         projection the one with the smallest tiltangle
%                         and set the 3D coordinate to 
%                         [markerset(1:2,irefmarkX,irefmark); imdim/2 + 1].
%   imdim               Dimension of images (default: 2048)
%  
%  OUTPUT
%   psi                 tilt axis azimuth in degrees.
%   tx, ty              N-vectors of the shift of the projections
%   sigma               RMS of linear fit
%   X                   The 3D coordinates of the markers. In contrast to
%                       TOM_ALIGNMENT3D the 3D coodinates are shifted, so
%                       that X(1:3,irefmark) == irefmarkX.
%
%See TOM_MARK_CVAF_ALIGNMENT3D_REPROJ how the output parameters determine
%the projection parameters, i.e. how the 3D points are projected onto their
%2D coordinates using psi, tx, ty and tiltangles.
%
%EXAMPLE
%   see the doTest function in the source code.
%
%REFERENCES
%
%SEE ALSO
%   TOM_ALIGNMENT3D, TOM_MARK_CVAF_ALIGNMENT3D_REPROJ
%
%   created by Thomas Haller, 09/14/2007
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




%if (0)
%    [psi, tx, ty, sigma, X]  = align3d(markerset, tiltangles, irefmark, irefmarkX, imdim);
%    doTest; return;
%else
%    [psi, tx, ty, sigma, X]  =  doTest_compare(markerset, tiltangles, irefmark, irefmarkX, imdim);
%end;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function [psi, tx, ty, sigma, X]  = align3d(markerset, tiltangles, irefmark, irefmarkX, imdim)


if (~exist('imdim','var') || isempty(imdim))
    imdim = 2048;
end;
imdim2 = imdim/2 + 1;


%   tilt nearest to zero
if (~exist('irefmark','var') || isempty(irefmark)) 
    irefmark = 1;
end;
if (~exist('irefmarkX','var') || isempty(irefmarkX))
    irefmarkX = find(abs(tiltangles) == min(abs(tiltangles)));
end;
if (length(irefmarkX) == 3)
elseif (length(irefmarkX)==1 && irefmarkX == round(irefmarkX))
    irefmarkX = [markerset(1:2, irefmarkX, irefmark)', imdim2];
else
    error('irefmarkX has an unexpected value.');
end;


ni = size(markerset, 2);
nm = size(markerset, 3);

markerset_defmask = all(isfinite(markerset), 1);



%   calculate means of difference vectors
meanx = zeros(nm, 1);
meany = zeros(nm, 1);
norm = zeros(nm,1);
for (imark = 1:nm)
    for (itilt = 1:ni)
        if (markerset_defmask(1,itilt,imark) && markerset_defmask(1,itilt,irefmark)) %allow overlapping MPs
            meanx(imark) = markerset(1,itilt,imark) - markerset(1,itilt,irefmark) + meanx(imark);
            meany(imark) = markerset(2,itilt,imark) - markerset(2,itilt,irefmark) + meany(imark);
            norm(imark) = norm(imark) +1;
        end;
    end;
end;
meanx = meanx ./ norm;
meany = meany ./ norm;




%   calculate some sums for determination of tilt axis azimuth
%   e.g. sumxx = sum(imark) sum(itilt) (delta(imark,itilt)-deltaaverage(imark))^2 
sumxx = 0;
sumyy = 0;
sumxy = 0;
for (imark = 1:nm)
    for (itilt = 1:ni)
        if (markerset_defmask(1,itilt,imark) && markerset_defmask(1,itilt,irefmark)) %allow overlapping MPs
            sumxx = (markerset(1,itilt,imark) - markerset(1,itilt,irefmark) - meanx(imark))^2 + sumxx;
            sumyy = (markerset(2,itilt,imark) - markerset(2,itilt,irefmark) - meany(imark))^2 + sumyy;
            sumxy = (markerset(1,itilt,imark) - markerset(1,itilt,irefmark) - meanx(imark)) * ...
                    (markerset(2,itilt,imark) - markerset(2,itilt,irefmark) - meany(imark)) + sumxy;
        end;
    end;
end;

%  calculate azimuth :
%  minimize sum( (x(i)*sin(psi) + y(i)*cos(psi))^2) =: Min(F)  --- linear regression  
%  d/d(psi) (F) leads to:
psi = 0.5*atan(2*sumxy/(sumxx-sumyy));
if (sumxx > sumyy)
    psi = psi - 0.5*pi*sign(psi);
end;



% calculate deviations
cpsi = cos(psi);
spsi = sin(psi);
%tpsi = tan(psi);
%sumt=0;
sumxx=0;
ndif =0;
Matrixmark4 = nan(1,ni,nm);
for (imark = 1:nm)
    for (itilt = 1:ni)
        if (markerset_defmask(1,itilt,imark) && markerset_defmask(1,itilt,irefmark)) %allow overlapping MPs
            if (imark ~= irefmark) 
                ndif = ndif +1; % count all markers except for refmark
            end;  
            Matrixmark4(1, itilt, imark) = (markerset(1,itilt,imark) - markerset(1,itilt,irefmark) - meanx(imark))*cpsi ...
                                         + (markerset(2,itilt,imark) - markerset(2,itilt,irefmark) - meany(imark))*spsi;
            sumxx = Matrixmark4(1,itilt,imark)^2 + sumxx;
        end;
    end;
end;
sigma = sqrt(sumxx/(ndif - nm));
%   deviation as angle in deg
%   sumtmp = fact*sqrt(sumtmp/ndif)*180/pi;
%disp(['Number of tilts:.............. = ' num2str(ni) ])
%disp(['minimum tilt: projection number= ' num2str(imintilt)]);
%disp(['Total number of marker points  = ' num2str(nm)])
%disp(['Number of reference point:.... = ' num2str(irefmark) ])
%disp(['Tilt axis azimuth:............ = ' num2str(psi*180/pi) ' deg'])
%disp(['RMS error fit:................ = ' num2str(sigma) '  pix'])





%   ---- 2nd part: determination of shifts ----

%   theta = tiltangle
theta_inrad = tiltangles * pi / 180; 
stheta = sin(theta_inrad);
ctheta = cos(theta_inrad);


%   determine 3D coordinates of markers -> solve eqn system

% x = zeros(size(Matrixmark,3)); y = zeros(size(Matrixmark,3)); z = zeros(size(Matrixmark,3)); 

%disp(['Coordinates of reference marker point ' num2str(irefmark) ' : ' ...
%        num2str(irefmarkX(1)) '  ' num2str(irefmarkX(2)) '  ' num2str(irefmarkX(3)) ]);
%disp(['Difference vectors of marker points:' ]);
X = nan(3,nm);
for (imark = 1:nm)
    sumxx = 0; 
    sumyy = 0; 
    sumxy = 0; 
    sumyx = 0; 
    salpsq = 0; 
    scalph = 0;
    P = zeros(3,3);
    temp = zeros(3);
    norm(:) = 0;
    for (itilt = 1:ni)
        if (markerset_defmask(1,itilt,imark) && markerset_defmask(1,itilt,irefmark)) %allow overlapping MPs
            norm(imark) = norm(imark) + 1;
            salpsq = salpsq + stheta(itilt)^2; %sum sin^2
            scalph = scalph + ctheta(itilt)*stheta(itilt); %sum cos*sin
            %sum delta x * cos
            sumxx = sumxx + (markerset(1,itilt,imark)-markerset(1,itilt,irefmark))*ctheta(itilt);
            %   sum delta(y)*cos
            sumyy = sumyy + (markerset(2,itilt,imark)-markerset(2,itilt,irefmark))*ctheta(itilt);
            %   sum delta(x)*sin
            sumxy = sumxy + (markerset(1,itilt,imark)-markerset(1,itilt,irefmark))*stheta(itilt);
            %   sum delta(y)*sin
            sumyx = sumyx + (markerset(2,itilt,imark)-markerset(2,itilt,irefmark))*stheta(itilt);
        end;
    end;
    P(1,1) = norm(imark) - salpsq*spsi^2;
    P(1,2) = salpsq*cpsi*spsi;
    P(1,3) = scalph*spsi;
    P(2,1) = P(1,2);
    P(2,2) = norm(imark) - salpsq*cpsi^2;
    P(2,3) = -scalph*cpsi;
    P(3,1) = P(1,3);
    P(3,2) = P(2,3);
    P(3,3) = salpsq;
    dt = det(P);
    temp(1) =   (sumxx*spsi-sumyy*cpsi)*spsi + ...
                (cpsi*meanx(imark)+spsi*meany(imark))*cpsi*norm(imark);
    temp(2) = - (sumxx*spsi-sumyy*cpsi)*cpsi ...
              + (cpsi*meanx(imark)+spsi*meany(imark)) * spsi*norm(imark);
    temp(3) = sumxy*spsi - sumyx*cpsi;
    if (dt ~= 0)
        P_t = P;
        P_t(1,1) = temp(1);
        P_t(2,1) = temp(2);
        P_t(3,1) = temp(3);
        X(1,imark) = det(P_t)/dt;
        P_t = P;
        P_t(1,2)=temp(1);
        P_t(2,2)=temp(2);
        P_t(3,2)=temp(3);
        X(2,imark) = det(P_t)/dt;
        P_t = P;
        P_t(1,3)=temp(1);
        P_t(2,3)=temp(2);
        P_t(3,3)=temp(3);
        X(3,imark) = det(P_t)/dt;
        %disp(['     ' num2str(imark) ' - ' num2str(irefmark) ' :.............. = ' ...
        %        num2str(X(1,imark)) '  ' num2str(X(2,imark)) '  ' num2str(X(3,imark)) ]);
        X(1,imark) = X(1,imark) + irefmarkX(1) - imdim/2 - 1;  % move to center
        X(2,imark) = X(2,imark) + irefmarkX(2) - imdim/2 - 1;
        X(3,imark) = X(3,imark) + irefmarkX(3) - imdim/2 - 1;
        %X(1,imark) = X(1,imark) + irefmarkX(1);  % move to center
        %X(2,imark) = X(2,imark) + irefmarkX(2);
        %X(3,imark) = X(3,imark) + irefmarkX(3);
    else
        %disp(['marker ' num2str(imark) ' : undefined! det = 0! Click more!']);
        markerset_defmask(1,:,imark) = false;
    end; 
end;




%   determination of shifts
Matrixmark56 = zeros(2, ni, nm);
Matrixmark910 = nan(2, ni, 1);
tx = nan(1, ni);
ty = nan(1, ni);
for (itilt = 1:ni)
    sumxx = 0;
    sumyy = 0;
    ndif = 0;
    for (imark = 1:nm)
        if (markerset_defmask(1,itilt,imark))
            shift1 = X(1,imark)*(spsi^2*ctheta(itilt)+cpsi^2) + ...
                     X(2,imark)*spsi*cpsi*(1-ctheta(itilt)) + ...
                     X(3,imark)*spsi*stheta(itilt) + imdim2;

            shift2 = X(1,imark)*spsi*cpsi*(1-ctheta(itilt)) + ...
                     X(2,imark)*(cpsi^2*ctheta(itilt)+spsi^2) ... 
                   - X(3,imark)*cpsi*stheta(itilt) + imdim2;
            ndif = ndif + 1;
            Matrixmark56(1,itilt,imark) = markerset(1,itilt,imark) - shift1;
            Matrixmark56(2,itilt,imark) = markerset(2,itilt,imark) - shift2;
            sumxx = sumxx + Matrixmark56(1,itilt,imark);
            sumyy = sumyy + Matrixmark56(2,itilt,imark);
        end;
    end;
    % mean values of shifts
    if (ndif > 0)
        tx(itilt) = sumxx / ndif;
        ty(itilt) = sumyy / ndif;
    else
        %Matrixmark(7:8,itilt,1) = 1000000;
    end;
    
    % deviations of individual shift from mean
    sumxx = 0;
    sumyy = 0;
    for (imark = 1:nm)
        if (markerset_defmask(1,itilt,imark))
            sumxx = sumxx + (Matrixmark56(1,itilt,imark)-tx(itilt))^2;
            sumyy = sumyy + (Matrixmark56(2,itilt,imark)-ty(itilt))^2;
        end;
    end;
    if (ndif > 1)
        Matrixmark910(1,itilt,1) = sqrt(sumxx/(ndif-1));
        Matrixmark910(2,itilt,1) = sqrt(sumyy/(ndif-1));
    else
        %Matrixmark(9:10,itilt,1)= -1000;
    end;
end;


X = X + imdim/2 + 1;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [psi, tx, ty, sigma, X]  = doTest_compare(markerset, tiltangles, irefmark, irefmarkX, imdim)


% Check whether the local function returns the same result as the
% original tom_alignment3d!

% tom_alignment3d ignores coodinates outside the image-dimension...
% thus, delete all those from the current markerset.
markerset_nan = markerset;
markerset_nan(:,~all(markerset>=1 & markerset<=imdim, 1)) = nan;
markerset_m1 = markerset_nan;
markerset_m1(~isfinite(markerset_nan)) = -1;

Matrixmark(1,1:length(tiltangles),1) = tiltangles;
Matrixmark(2:3,1:size(markerset_m1,2), 1:size(markerset_m1,3)) = markerset_m1;

nrepeat = 1;
tic;
for (i=1:nrepeat)
    [psi, tx, ty, sigma, X]  = align3d(markerset_nan, tiltangles, irefmark, irefmarkX, imdim);
end;
t1 = toc;
irefmarkX = X(:,irefmark);
tic;
for (i=1:nrepeat)
    %[Matrixmark_2, psi_2, sigma_2, x_2, y_2, z_2]  = tom_alignment3d(Matrixmark, irefmark, 0, X(:,irefmark), imdim);
    [Matrixmark_2, psi_2, sigma_2, x_2, y_2, z_2]  = tom_Rec3dRigidBodyAlignment(Matrixmark, irefmark, 0, irefmarkX, imdim, 'off');
end;
t2 = toc;

X_2 = [x_2;y_2;z_2] + (imdim/2+1);
disp(['original: ' num2str(t2) ' sek; new: ' num2str(t1)]);
maxerror = 1e3*eps;
if (max(max(abs(X - X_2))) > maxerror || ...
    abs(psi - psi_2) > maxerror || ...
    abs(sigma_2 - sigma) > maxerror || ...
    max(max(abs((Matrixmark_2(7:8,:,1) - [tx;ty])))) > maxerror)
    warning('Different results!!!!!');
end;        



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function doTest

imdim = 1024;

rand('state', 3);
if (0)
    Matrixmark = tom_emread('/fs/bmsan/apps/tom_dev/Reconstruction/tom_mark/test_data/dualaxis/DUAL/marker4_pyMOD.em');
    tiltangles = Matrixmark.Value(1,:,1);
    ni = size(Matrixmark.Value, 2);
    nm = size(Matrixmark.Value, 3);
    
    markerset = Matrixmark.Value(2:3,:,:);
    Matrixmark = Matrixmark.Value;

    idx_undef = ~all(markerset >= 1 & markerset <= imdim, 1);
    
    imintilt = find(abs(tiltangles) == min(abs(tiltangles)), 1);
    irefmark = find(~idx_undef(1,imintilt,:), 1);
    
    irefmark3D = [markerset(:, imintilt, irefmark); imdim/2+1];
    psi = nan;
    txy = nan(2,ni);
    X = nan(3,nm);
    
else
    ni = 100;
    nm = 4;

    tiltangles = linspace(-60,60,ni);
    imintilt = find(abs(tiltangles) == min(abs(tiltangles)), 1);
    
    psi = rand() * 2*pi;
    if (psi > pi)
        psi = psi - 2*pi;
    end;
    txy = randn(2, ni) * imdim/50;
    txy(1:2, imintilt) = 0;

    X = rand(3, nm) .* repmat([imdim,imdim,imdim/2]',[1,nm]);

    [P, markerset] = tom_mark_cvaf_alignment3d_reproj(tiltangles, psi, txy(1,:), txy(2,:), imdim, X);

    idx_undef = ~all(markerset >= 1 & markerset <= imdim, 1);
    irefmark = find(~idx_undef(1,imintilt,:), 1);
    
    Matrixmark(1,:,1) = tiltangles;
    Matrixmark(2:3,1:ni,1:nm) = markerset;

    irefmark3D = X(1:3,irefmark);

end;

markerset(1:2, idx_undef) = nan;
Matrixmark(2:3, idx_undef) = -1;




[alignment3d_Matrixmark, alignment3d_psi, alignment3d_sigma, alignment3d_x, alignment3d_y, alignment3d_z]  = tom_alignment3d(Matrixmark, irefmark, imintilt, irefmark3D, imdim);
alignment3d_X = [alignment3d_x; alignment3d_y; alignment3d_z];



[align3d_psi, align3d_tx, align3d_ty, align3d_sigma, align3d_X]  = align3d(markerset, tiltangles, irefmark, irefmark3D, imdim);


alignment3d_X = alignment3d_X + (imdim/2+1);
align3d_X = align3d_X + (imdim/2+1);

%TODO: check tx, ty

disp(['psi = ' num2str(psi) ' ; ' num2str(alignment3d_psi) ' ; ' num2str(align3d_psi)]);
disp(['sigma = ----- ; ' num2str(alignment3d_sigma) ' ; ' num2str(align3d_sigma)]);

disp(['max dist mean shift: orig vs. tom_alignment3d:                      ' num2str(max(sqrt(sum((alignment3d_Matrixmark(7:8,:,1) - [txy]) .^ 2,1))))])
disp(['max dist mean shift: orig vs. tom_mark_cvaf_alignment3d:            ' num2str(max(sqrt(sum(([align3d_tx;align3d_ty] - [txy]) .^ 2,1))))])
disp(['max dist mean shift: tom_mark_cvaf_alignment3d vs. tom_alignment3d: ' num2str(max(sqrt(sum((alignment3d_Matrixmark(7:8,:,1) - [align3d_tx;align3d_ty]) .^ 2,1))))])
disp(['max dist X: orig vs. tom_align3d:                                   ' num2str(max(sqrt(sum((X - align3d_X).^2,1))))]);
disp(['max dist X: orig vs. tom_alignment3d:                               ' num2str(max(sqrt(sum((alignment3d_X - X).^2,1))))]);
disp(['max dist X: tom_align3d vs. tom_alignment3d:                        ' num2str(max(sqrt(sum((alignment3d_X - align3d_X).^2,1))))]);

