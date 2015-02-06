function integ_ps = tom_image2ctf(image,logscale,centerregion,split,nointeg)
%TOM_IMAGE2CTF creates ...
%
%   integ_ps = tom_image2ctf(image,logscale,centerregion,split)
%
%PARAMETERS
%
%  INPUT
%   image               ...
%   logscale            ...
%   centerregion        ...
%   split               ...
%  
%  OUTPUT
%   data		...
%
%EXAMPLE
%   .... = tom_image2ctf(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by ... (author date)
%   updated by ...
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

error(nargchk(0, 5, nargin, 'struct'))

if nargin < 5
    nointeg = 0;
end

if nargin < 4
    split = 1;
end

if nargin < 3
    centerregion = 0;
end

if nargin < 2
    logscale = 1;
end

if split > 1
    ps_out = zeros(size(image)./split);
    for i=1:split
        if (i==1) 
            im_old=image;
        end
        image = split_image(im_old,split,i);
        image = tom_smooth(image,round(size(image,1).*.05)); % 5 % border smoothing
        ps = tom_ps(image); 
        ps_out = ps_out + ps;    
    end
else
    ps_out = tom_ps(image);
end

ps = ps_out;

if centerregion == 1
    %cut out central part of power spectrum
    size_ps = size(ps);
    center = floor(size_ps ./ 2 + 1);
    width = floor(center ./ 2 - 1);
    ps_out = ps(center(1)-width(1):center(1)+width(1)-1,center(2)-width(2):center(2)+width(2)-1);
else
    size_ps = size(ps);
    center = floor(size_ps ./ 2 + 1);
    width = floor(center - 1);
    ps_out = ps;
end

%circular integration over ps
if nointeg ~= 1
    integ_ps = tom_cart2polar(ps_out);
    integ_ps = sum(integ_ps,2)./(size(integ_ps,2)); % fix me ???? or size(integ_ps,1)
else
    integ_ps = ps_out;
end

if logscale == 1
    integ_ps = log(integ_ps);
end

%integ_ps = integ_ps';
if nointeg ~= 1
    %integ_ps = integ_ps';
        integ_ps = smooth(integ_ps,round(length(integ_ps).*.05))';
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  split image for ps averaging                                       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function split_image=split_image(im,number_of_splits,split_nr)
im_sz=size(im,1);

inkre=round(im_sz./number_of_splits);
start=((split_nr-1).*inkre)+1;
stop=((split_nr).*inkre);
split_image=im(start:stop,start:stop);
