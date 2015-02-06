function [Matrixmark, psi, sigma, x, y, z]  = tom_alignment3dDisp ...
(Matrixmark,irefmark,imintilt,r,imdim,DisplayOnOff)
%TOM_ALIGNMENT3D calculates the tilt-axis, shifts, fit-errors of an EM markerfile
%
%   [Matrixmark, psi, sigma, x, y, z]  = tom_alignment3dDisp(Matrixmark,irefmark,imintilt,r,imdim,DisplayOnOff)
%
% TOM_ALIGNMENT3D calculates the tilt-axis, the shifts and fit-errors of an
% EM markerfile
%
%PARAMETERS
%
%  INPUT
%   Matrixmark[]        Matrix of dimension (10, no. of images, no. of markers)
%   irefmark            Reference marker (default: marker 1)
%   imintilt            reference projection (default: Projection with minimum tilt)
%   r[3]                3D coordinates assigned to marker irefmark (default: r = [x(imintilt), y(imintilt), imdim/2 + 1])
%   imdim               Dimension of images (default: 2048)
%  
%  OUTPUT
%   Matrixmark          as input but altered:
%                       Matrixmark contains 10 columns:
%                       1          2          3          4          5          6          7          8          9          10
%                     tilt       x-coord    y-coord   dev of mark shift(x)  shift(y)    mean-shift mean-shift dev shift  dev shift   
%
%                       additionally : Matrixmark(7,1,2)=psi (in deg);
%
%   psi                 tilt axis azimuth       
%   sigma               RMS of linear fit
%   x[no. of markers]   3D-coordinates of markers
%   y[no. of markers]
%   z[no. of markers]
%
%EXAMPLE
%   ... = tom_alignment3dDisp(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   tom_marker, tom_reconstruction3d
%
%   created by 07/25/02
%   updated by FF 08/02/03
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

%   tilt nearest to zero
error(nargchk(1,6,nargin))
if (nargin < 3)
    [dummy, imintilt]=min(abs(Matrixmark(1,:,1)));
end;
%  
if (nargin < 2)
    irefmark = 1;
end;
if (nargin < 5)
    imdim = 2048;
end;
if (nargin < 4)
    %  define r(1), r(2), and r(3) of irefmark
    r = [Matrixmark(2,imintilt,irefmark) Matrixmark(3,imintilt,irefmark) (imdim/2 +1)];
end;
if (nargin < 6)
    DisplayOnOff = 'off';
end;
    

%   calculate means of difference vectors
meanx=zeros(size(Matrixmark,3),1);
meany=zeros(size(Matrixmark,3),1);
norm=zeros(size(Matrixmark,3),1);
for imark = 1:size(Matrixmark,3),
    for itilt = 1:size(Matrixmark,2),
        if (Matrixmark(2,itilt,imark)> -1.0) & (Matrixmark(2,itilt,irefmark)> -1.0)%allow overlapping MPs
            meanx(imark) = Matrixmark(2,itilt,imark) - Matrixmark(2,itilt,irefmark) + meanx(imark);
            meany(imark) = Matrixmark(3,itilt,imark) - Matrixmark(3,itilt,irefmark) + meany(imark);
            norm(imark) = norm(imark) +1;
        end;
    end;
end;
meanx(:) = meanx(:) ./ norm(:);
meany(:) = meany(:) ./ norm(:);
%   calculate some sums for determination of tilt axis azimuth
%   e.g. sumxx = sum(imark) sum(itilt) (delta(imark,itilt)-deltaaverage(imark))^2 


sumxx=0;
sumyy=0;
sumxy=0;
sumtmp=0;
for imark = 1:size(Matrixmark,3),
    for itilt = 1:size(Matrixmark,2),
        if (Matrixmark(2,itilt,imark)> -1.0) & (Matrixmark(2,itilt,irefmark)> -1.0)%allow overlapping MPs
            sumxx = (Matrixmark(2,itilt,imark) - Matrixmark(2,itilt,irefmark)-meanx(imark))^2 + sumxx;
            sumyy = (Matrixmark(3,itilt,imark) - Matrixmark(3,itilt,irefmark)-meany(imark))^2 + sumyy;
            sumxy = (Matrixmark(2,itilt,imark) - Matrixmark(2,itilt,irefmark)-meanx(imark)) * (Matrixmark(3,itilt,imark) - Matrixmark(3,itilt,irefmark)-meany(imark)) + sumxy;
%            sumtmp  = (Matrixmark(2,itilt,imark) - Matrixmark(2,itilt,irefmark))^2;
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
psiindeg = psi*180/pi;

%  psi 2?!
fact = 1 /(1 + (sumxy/sumxx)^2);

% calculate deviations
cpsi = cos(psi);
spsi = sin(psi);
tpsi = tan(psi);
sumt=0;
sumxx=0;
ndif =0;
for imark = 1:size(Matrixmark,3),
    for itilt = 1:size(Matrixmark,2),
        if (Matrixmark(2,itilt,imark)> -1.0) & (Matrixmark(2,itilt,irefmark)> -1.0)%allow overlapping MPs
            if (imark ~= irefmark) 
                ndif= ndif +1; % count all markers except for refmark
            end;  
            Matrixmark(4,itilt,imark)= (Matrixmark(2,itilt,imark) - Matrixmark(2,itilt,irefmark)-meanx(imark))*cpsi ...
                + (Matrixmark(3,itilt,imark) - Matrixmark(3,itilt,irefmark)-meany(imark))*spsi;
            sumxx=Matrixmark(4,itilt,imark)^2 +sumxx;
        end;
    end;
end;
sigma = sqrt(sumxx/(ndif - size(Matrixmark,3)));
%   deviation as angle in deg
%   sumtmp = fact*sqrt(sumtmp/ndif)*180/pi;

if (strcmp(DisplayOnOff,'on')) == 1
disp(['Number of tilts:.............. = ' num2str(size(Matrixmark,2)) ])
disp(['minimum tilt: projection number= ' num2str(imintilt)]);
disp(['Total number of marker points  = ' num2str(size(Matrixmark,3))])
disp(['Number of reference point:.... = ' num2str(irefmark) ])
disp(['Tilt axis azimuth:............ = ' num2str(psiindeg) ' deg'])
disp(['RMS error fit:................ = ' num2str(sigma) '  pix'])
end;

%   ---- 2nd part: determination of shifts ----

%   theta = tiltangle
theta_inrad=squeeze(squeeze(Matrixmark(1,:,1)))*pi/180; 
stheta = sin(theta_inrad(:)) ;
ctheta = cos(theta_inrad(:)) ;
%Matrixmark(5,:,1) = stheta;
%Matrixmark(6,:,1) = ctheta;


%   determine 3D coordinates of markers -> solve eqn system

% x = zeros(size(Matrixmark,3)); y = zeros(size(Matrixmark,3)); z = zeros(size(Matrixmark,3)); 

if (strcmp(DisplayOnOff,'on')) == 1
disp(['Coordinates of reference marker point ' num2str(irefmark) ' : ' ...
        num2str(r(1)) '  ' num2str(r(2)) '  ' num2str(r(3)) ]);
disp(['Difference vectors of marker points:' ]);
end;

for imark = 1:size(Matrixmark,3),
    sumxx=0; sumyy=0; sumxy=0; sumyx=0; salpsq = 0; scalph = 0;
    P = zeros(3,3);
    temp = zeros(3);
    norm(:)=0;
    for itilt = 1:size(Matrixmark,2),
        if (Matrixmark(2,itilt,imark)> -1.0) & (Matrixmark(2,itilt,irefmark)> -1.0)%allow overlapping MPs
            norm(imark)= norm(imark)+1;
            salpsq = salpsq + stheta(itilt)^2; %sum sin^2
            scalph = scalph + ctheta(itilt)*stheta(itilt); %sum cos*sin
            %sum delta x * cos
            sumxx = sumxx + (Matrixmark(2,itilt,imark) - Matrixmark(2,itilt,irefmark))*ctheta(itilt);
            %   sum delta(y)*cos
            sumyy = sumyy + (Matrixmark(3,itilt,imark) - Matrixmark(3,itilt,irefmark))*ctheta(itilt);
            %   sum delta(x)*sin
            sumxy = sumxy + (Matrixmark(2,itilt,imark) - Matrixmark(2,itilt,irefmark))*stheta(itilt);
            %   sum delta(y)*sin
            sumyx = sumyx + (Matrixmark(3,itilt,imark) - Matrixmark(3,itilt,irefmark))*stheta(itilt);
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
    temp(1) = (sumxx*spsi-sumyy*cpsi)*spsi + ...
        (cpsi*meanx(imark)+spsi*meany(imark))*cpsi*norm(imark);
    temp(2) = -(sumxx*spsi-sumyy*cpsi)*cpsi + (cpsi*meanx(imark)+spsi*meany(imark))*...
        spsi*norm(imark);
    temp(3) = sumxy*spsi - sumyx*cpsi;
    if (dt ~= 0)
        P_t=P;
        P_t(1,1) = temp(1);
        P_t(2,1) = temp(2);
        P_t(3,1) = temp(3);
        x(imark) = det(P_t)/dt;
        P_t = P;
        P_t(1,2)=temp(1);
        P_t(2,2)=temp(2);
        P_t(3,2)=temp(3);
        y(imark) = det(P_t)/dt;
        P_t = P;
        P_t(1,3)=temp(1);
        P_t(2,3)=temp(2);
        P_t(3,3)=temp(3);
        z(imark) = det(P_t)/dt;
        
        if (strcmp(DisplayOnOff,'on')) == 1
        disp(['     ' num2str(imark) ' - ' num2str(irefmark) ' :.............. = ' ...
                num2str(x(imark)) '  ' num2str(y(imark)) '  ' num2str(z(imark)) ]);
        end;
        x(imark) = x(imark) + r(1) - imdim/2 -1;  % move to center
        y(imark) = y(imark) + r(2) - imdim/2 -1;
        z(imark) = z(imark) + r(3) - imdim/2 -1;
    else

        if (strcmp(DisplayOnOff,'on')) == 1
        disp(['marker ' num2str(imark) ' : undefined! det = 0! Click more!']);
        end;

        x(imark) = 1000000;
    end; 
end;

%   determination of shifts

shift=zeros(2,size(Matrixmark,2));
for itilt = 1:size(Matrixmark,2),
    sumxx = 0;
    sumyy = 0;
    ndif = 0;
    for imark = 1:size(Matrixmark,3),
        shift(1,itilt) = x(imark)*(spsi^2*ctheta(itilt)+cpsi^2) + ...
                         y(imark)*spsi*cpsi*(1-ctheta(itilt)) + ...
                         z(imark)*spsi*stheta(itilt) + imdim/2 + 1;

        shift(2,itilt) = x(imark)*spsi*cpsi*(1-ctheta(itilt)) + ...
                         y(imark)*(cpsi^2*ctheta(itilt)+spsi^2) ... 
                       - z(imark)*cpsi*stheta(itilt) + imdim/2 + 1;
        if (Matrixmark(2,itilt,imark)> -1.0) & (x(imark)~=1000000)
            ndif = ndif + 1;
            Matrixmark(5,itilt,imark) = Matrixmark(2,itilt,imark) - shift(1,itilt);
            Matrixmark(6,itilt,imark) = Matrixmark(3,itilt,imark) - shift(2,itilt);
           % shift(1,itilt) = Matrixmark(5,itilt,imark);
           % shift(2,itilt) = Matrixmark(6,itilt,imark);
        else
            Matrixmark(5:6,itilt,imark) = 0;
            %shift(:,itilt) = 0;
        end;
        sumxx = sumxx + Matrixmark(5,itilt,imark);
        sumyy = sumyy + Matrixmark(6,itilt,imark);
    end;
    % mean values of shifts
    if (ndif > 0)
        Matrixmark(7,itilt,1) = sumxx / ndif;
        Matrixmark(8,itilt,1) = sumyy / ndif;
    else
        Matrixmark(7:8,itilt,1) = 1000000;
    end;
    
    % deviations of individual shift from mean
    sumxx = 0;
    sumyy = 0;
    for imark = 1:size(Matrixmark,3),
        if (Matrixmark(2,itilt,imark)> -1.0) & (x(imark)~=1000000)
            sumxx = sumxx + (Matrixmark(5,itilt,imark)-Matrixmark(7,itilt,1))^2;
            sumyy = sumyy + (Matrixmark(6,itilt,imark)-Matrixmark(8,itilt,1))^2;
        end;
    end;
    if (ndif > 1)
        Matrixmark(9,itilt,1)=sqrt(sumxx/(ndif-1));
        Matrixmark(10,itilt,1)=sqrt(sumyy/(ndif-1));
    else
        Matrixmark(9:10,itilt,1)= -1000;
    end;
end;
% store tilt axis somewhere (delete?!)
Matrixmark(7,1,2)=psiindeg;