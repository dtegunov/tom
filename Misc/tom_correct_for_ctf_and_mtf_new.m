function [out ctf_theory ctf_correct ctf_flip]=tom_correct_for_ctf_and_mtf_new(in,Fit,method,cutoff,mtf)

%tom_correct_for_ctf_and_mtf_new  Corrects a volume for ctf and mtf
%
%   out=tom_correct_for_ctf_and_mtf_new(in,Fit,method,cutoff)
%
%
%PARAMETERS
%
%  INPUT
%   in                  original image
%   Fit                 Fit structure from tom_fit_ctf_gui.m
%   method              method for phase flipping etc.
%   cutoff              cutoff frequency in pixel
%   mtf                 MTF of the CCD camera.
%
%  OUTPUT
%   corrected           corrected volume
%
%EXAMPLE
%       out=tom_correct_for_ctf_and_mtf_new(in,Fit,'phase&mtf',32,mtf)
%
%
%REFERENCES
%
%SEE ALSO
%
%   created by SN 01/07/09
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

if ~isstruct(Fit)
    error('Fit structure from tom_fit_ctf_gui.m is required.');
end

if nargin < 4
    cutoff = size(in,1)./2
end

if nargin < 3
    method = 'flip';
end

out=0;


Dz=Fit.Dz_det;
Dz_delta=Fit.Dz_delta_det;
Phi_0=Fit.Phi_0_det;
pix_size=Fit.EM.Objectpixelsize;
voltage=Fit.EM.Voltage;
img_size=size(in);
sigma=Fit.decay_det;
q_0_r=Fit.decay_part_coh_ill_det;
Cs=Fit.EM.Cs;
Cc=Fit.EM.Cc;
deltaE=Fit.decay_energy_spread_det;
amplitude_contrast=Fit.amplitude_contrast_det;

[phase amplitude decay E_k E_i_k E_e_k]=tom_ctf2(Dz,Dz_delta,Phi_0,pix_size,voltage,img_size,Cs,sigma,q_0_r,Cc,deltaE);

if Fit.decay_det==0
    decay=1;
end;

ctf_theory=E_e_k.*E_i_k.*decay.*(sqrt(1-amplitude_contrast.^2).*phase + amplitude_contrast.*amplitude);

ctf_correct=tom_limit(ctf_theory./abs(ctf_theory+0.00001),-1,1);

mask = tom_spheremask(ones(size(ctf_correct)),cutoff,cutoff./32);

if strcmp(method,'flip&mtf')
    disp('flip&mtf');
    mtf_img = tom_polar2cart(tom_norm(mtf,1));
    mtf_img = imresize(mtf_img,[size(ctf_correct,1) size(ctf_correct,2)],'bicubic');
    ctf_correct=ctf_correct.*(1./(sqrt(mtf_img)+abs(1-mask)));
end;

if strcmp(method,'flip')
    %disp('flip');
    %mtf_img = tom_polar2cart(tom_norm(mtf,1));
    %mtf_img = imresize(mtf_img,[size(ctf_correct,1) size(ctf_correct,2)],'bicubic');
    %ctf_correct=ctf_correct.*(1./(mtf_img+abs(1-mask)));
    ctf_correct=ctf_correct;
end;

if strcmp(method,'wiener&mtf')
%    disp('wiener&mtf');
    if size(mtf,2)>1
        mtf_img = mtf;
        if size(mtf_img)~=size(ctf_correct)
            mtf_img = imresize(mtf_img,[size(ctf_correct,1) size(ctf_correct,2)],'bicubic');
        end;        
    else
        mtf_img = tom_polar2cart(tom_norm(mtf,1));
        mtf_img = imresize(mtf_img,[size(ctf_correct,1) size(ctf_correct,2)],'bicubic');
        ctf_correct=ctf_correct.*(1./(sqrt(mtf_img)+abs(1-mask)));
    end;
   
    
    
end;


ctf_correct=ctf_correct.*mask;

%in=tom_smooth(in,32,10,'zero');
ctf_flip=ctf_theory.*ctf_correct;

if strcmp(method,'wiener&mtf')
    
    out = deconvwnr(in,double(real(fftshift(ifft2(fftshift(ctf_theory.*(1./sqrt(mtf_img))))))),2);
    out = tom_apply_weight_function(out,mask);
else
    out = tom_apply_weight_function(in,ctf_correct,mask);
end;

