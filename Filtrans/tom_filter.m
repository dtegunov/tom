function im = tom_filter(im,radius,flag,space,lambda)
%TOM_FILTER convolutes with spherical or quadratic kernel
%
%   im = tom_filter(im,radius,flag,space,lambda)
%
%PARAMETERS
%
%  INPUT
%   im                  one, two or 3D volume
%   radius              radius of kernel for spherical kernel for quadratic kernel: side length of kernel
%   flag                'circ' for spherical kernel, 'quadr' for quadratic kernel
%   space               'real' for filtering in real space, 'fourier' for filtering in fourier space (default)
%   lambda              ...
%  
%  OUTPUT
%   im                  filtered image
%
%EXAMPLE
% Usage for anisotropic diffusion
% ===============================
%
%  filtered = tom_filter(im,function,sigma,iterations,lambda)
% 
%   function: influence function that determines how the
%      diffusion depends on the local image gradient.  Three
%      example psi functions are:
%         aniso_linearPsi
%         aniso_lorentzianPsi
%         aniso_tukeyPsi
%      Default 'aniso' is tukeyPsi, that makes this the same as one of
%      original two algorithms proposed by Perona et al.  But
%      tukeyPsi is often a better choice because it does a better
%      job of maintaining the sharpness of an edge.  LinearPsi
%      gives standard linear diffusion, i.e., shift-invariant
%      convolution by a Gaussian. 
%   sigma: scale parameter on the psiFunction.  Choose this
%      number to be bigger than the noise but small than the real
%      discontinuties. [Default = 1]
%   iterations: number of iterations [Default = 10]
%   lambda: rate parameter [Default = 1/4] To approximage a
%      continuous-time PDE, make lambda small and increase the
%      number of iterations.
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by AK 04/01/06 added real space filtering AK
%   updated by AK 06/03/06 added anisotropic filtering AK
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

error(nargchk(2,5,nargin));
if (nargin < 3)
    flag = 'circ';
end

if (nargin < 4)
    space = 'fourier';
end

%anisotropic diffusion
if strcmp(radius,'aniso') == 1 | strcmp(radius,'aniso_linearPsi') == 1 | strcmp(radius,'aniso_lorentzianPsi') == 1 | strcmp(radius,'aniso_tukeyPsi') == 1
    
    
    if strcmp(radius,'aniso_lorentzianPsi')
        fun = 'lorentzianPsi';
    elseif strcmp(radius,'aniso_linearPsi')
        fun = 'linearPsi';
    else
        fun = 'tukeyPsi';
    end

    if nargin < 3
        sigma = 1;
    else
        sigma = flag;
    end
    
    if nargin < 4
        iterations = 10;
    else 
        iterations = space;
    end
    
    if nargin < 5
        lambda = 0.25;
    end
    
    if size(im,3) == 1
        im = aniso(double(im),fun,sigma,iterations,lambda);
    else
        im = aniso3(double(im),fun,sigma,iterations,lambda);
    end
    
    
else   


    if strcmp(lower(space),'fourier') == 1
        %filtering in Fourier space

        switch lower(flag)
            case 'circ'
                mask=ones(size(im,1),size(im,2),size(im,3));
                mask=tom_spheremask(mask,radius);
            case 'quadr'
                mask=zeros(size(im,1),size(im,2),size(im,3));
                if size(im,3) > 1
                    cent=[floor(size(im,1)/2)+1 floor(size(im,2)/2)+1 floor(size(im,3)/2)+1];
                    mask(cent(1)-floor(radius/2):cent(1)-floor(radius/2)+radius-1,...
                        cent(2)-floor(radius/2):cent(2)-floor(radius/2)+radius-1,cent(3)-floor(radius/2):cent(3)-floor(radius/2)+radius-1)=1;
                elseif ((size(im,3) == 1) &  (size(im,2) > 1))
                    cent=[floor(size(im,1)/2)+1 floor(size(im,2)/2)+1];
                    mask(cent(1)-floor(radius/2):cent(1)-floor(radius/2)+radius-1,...
                        cent(2)-floor(radius/2):cent(2)-floor(radius/2)+radius-1)=1;
                else
                    cent=[floor(size(im,1)/2)+1];
                    mask(cent(1)-floor(radius/2):cent(1)+radius-1)=1;
                end
        end

        npix=sum(sum(sum(mask)));
        im=real(fftshift(tom_ifourier(tom_fourier(mask).*tom_fourier(im))/npix));

    else
        %filtering in Real Space

        switch lower(flag)
            case 'circ'
                if size(im,3) == 1
                    %2D
                    h = fspecial('disk',radius);
                else
                    %3D
                    error('Not implemented.');
                end
            case 'quadr'
                if size(im,3) == 1
                    %2D
                    h = fspecial('average',[radius radius]);
                else
                    %3D
                    h = ones([radius,radius,radius])./prod([radius,radius,radius]);
                end
        end

        im = imfilter(im,h);
    end
end