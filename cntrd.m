function out=cntrd(im,mx,sz,interactive)
% out=cntrd(im,mx,sz,interactive)
% 
% PURPOSE:  calculates the centroid of bright spots to sub-pixel accuracy.
%  Inspired by Grier & Crocker's feature for IDL, but greatly simplified and optimized
%  for matlab
% 
% INPUT:
% im: image to process, particle should be bright spots on dark background with little noise
%   ofen an bandpass filtered brightfield image or a nice fluorescent image
%
% mx: locations of local maxima to pixel-level accuracy from pkfnd.m
%
% sz: diamter of the window over which to average to calculate the centroid.  
%     should be big enough
%     to capture the whole particle but not so big that it captures others.  
%     if initial guess of center (from pkfnd) is far from the centroid, the
%     window will need to be larger than the particle size.  RECCOMMENDED
%     size is the long lengthscale used in bpass plus 2.
%     
%
% interactive:  OPTIONAL INPUT set this variable to one and it will show you the image used to calculate  
%    each centroid, the pixel-level peak and the centroid
%
% NOTE:
%  - if pkfnd, and cntrd return more then one location per particle then
%  you should try to filter your input more carefully.  If you still get
%  more than one peak for particle, use the optional sz parameter in pkfnd
%  - If you want sub-pixel accuracy, you need to have a lot of pixels in your window (sz>>1). 
%    To check for pixel bias, plot a histogram of the fractional parts of the resulting locations
%  - It is HIGHLY recommended to run in interactive mode to adjust the parameters before you
%    analyze a bunch of images.
%
% OUTPUT:  a N x 4 array containing, x, y and brightness for each feature
%           out(:,1) is the x-coordinates
%           out(:,2) is the y-coordinates
%           out(:,3) is the brightnesses
%           out(:,4) is the sqare of the radius of gyration
%
% CREATED: Eric R. Dufresne, Yale University, Feb 4 2005
%  5/2005 inputs diamter instead of radius
%  Modifications:
%  D.B. (6/05) Added code from imdist/dist to make this stand alone.
%  ERD (6/05) Increased frame of reject locations around edge to 1.5*sz
%  ERD 6/2005  By popular demand, 1. altered input to be formatted in x,y
%  space instead of row, column space  2. added forth column of output,
%  rg^2
%  ERD 8/05  Outputs had been shifted by [0.5,0.5] pixels.  No more!
%  ERD 8/24/05  Woops!  That last one was a red herring.  The real problem
%  is the "ringing" from the output of bpass.  I fixed bpass (see note),
%  and no longer need this kludge.  Also, made it quite nice if mx=[];
%  ERD 6/06  Added size and brightness output ot interactive mode.  Also 
%   fixed bug in calculation of rg^2
%  JWM 6/07  Small corrections to documentation 


if nargin==3
   interactive=0; 
end

if sz/2 == floor(sz/2)
warning('sz must be odd, like bpass');
end

if isempty(mx)
    warning('there were no positions inputted into cntrd. check your pkfnd theshold')
    out=[];
    return;
end


r=(sz+1)/2;
%create mask - window around trial location over which to calculate the centroid
m = 2*r;
x = 0:(m-1) ;
cent = (m-1)/2;
x2 = (x-cent).^2;
dst=zeros(m,m);
for i=1:m
    dst(i,:)=sqrt((i-1-cent)^2+x2);
end


ind=find(dst < r);

msk=zeros([2*r,2*r]);
msk(ind)=1.0;
%msk=circshift(msk,[-r,-r]);

dst2=msk.*(dst.^2);
ndst2=sum(sum(dst2));

[nr,nc]=size(im);
%remove all potential locations within distance sz from edges of image
ind=find(mx(:,2) > 1.5*sz & mx(:,2) < nr-1.5*sz);
mx=mx(ind,:);
ind=find(mx(:,1) > 1.5*sz & mx(:,1) < nc-1.5*sz);
mx=mx(ind,:);

[nmx,crap] = size(mx);

%inside of the window, assign an x and y coordinate for each pixel
xl=zeros(2*r,2*r);
for i=1:2*r
    xl(i,:)=(1:2*r);
end
yl=xl';

pts=[];
%loop through all of the candidate positions
for i=1:nmx
    %create a small working array around each candidate location, and apply the window function
    tmp=msk.*im((mx(i,2)-r+1:mx(i,2)+r),(mx(i,1)-r+1:mx(i,1)+r));
    %calculate the total brightness
    norm=sum(sum(tmp));
    %calculate the weigthed average x location
    xavg=sum(sum(tmp.*xl))./norm;
    %calculate the weighted average y location
    yavg=sum(sum(tmp.*yl))./norm;
    %calculate the radius of gyration^2
    %rg=(sum(sum(tmp.*dst2))/ndst2);
    rg=(sum(sum(tmp.*dst2))/norm);
    
    %concatenate it up
    pts=[pts,[mx(i,1)+xavg-r,mx(i,2)+yavg-r,norm,rg]'];
    
    %OPTIONAL plot things up if you're in interactive mode
    if interactive==1
     imagesc(tmp)
     axis image
     hold on;
     plot(xavg,yavg,'x')
     plot(xavg,yavg,'o')
     plot(r,r,'.')
     hold off
     title(['brightness ',num2str(norm),' size ',num2str(sqrt(rg))])
     pause
    end

    
end
out=pts';

