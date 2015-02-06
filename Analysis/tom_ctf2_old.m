function [phase amplitude decay E_k E_i_k E_e_k]=tom_ctf2(Dz,Dz_delta,Phi_0,pix_size,voltage,img_size,Cs,sigma,q_0_r,Cc,deltaE)
%TOM_CTF2D calculates 2D contrast transfer function (phase only!),
%           including astigmatism.
%
%   [phase amplitude decay E_k E_i_k E_e_k]==tom_ctf2(Dz,Dz_delta,Phi_0,pix_size,voltage,img_size,Cs)
%
%PARAMETERS
%  INPUT
%   Dz:         defocus of the objective lens in m (<0 underfocus, >0 overfocus).
%   Dz_delta:   the focal difference due to axial astigmatism in m.
%   Phi_0:      the reference angle defining the azimuthal direction of the
%               axial astigmatism in degree.
%   pix_size:   pixelsize in the object-plane in m.
%   voltage:    accelerating voltage in V.
%   img_size:   image size of newlz created image in pixel.
%   Cs:         the third-order spherical aberration constant in m.
%   sigma:      envelope of ctf (optional). If chosen decay 
%               ctf ~exp(-(freq/sigma)^2) is assumed. SIGMA is in units of Nyqvist.
%               => ctf(sigma) = 1/e 
%   q_0_r:      a quantity of dimension 1/length specifying the size of the
%               source in reduced coordinates, see Frank page 43, 2.19a
%               try 0.01-0.05
%   Cc:         the fourth-order chromatic aberration constant in m,
%               typically 2.2e-3 m.
%   deltaE:     energy spread, eg. 0.8;
%
%  OUTPUT
%   phase:      output containing the centro-symmetric phase. 
%   amplitude:  output containing the centro-symmetric amplitude. 
%   decay:      envelope function with an exponential decay.
%   E_k:        'compound envelope function', Frank formula 2.14: E_k=E_i_k * E_e_k
%   E_i_k:      envelope function based on partially coherent illumination
%   E_e_k:      envelope function due to energy spread 
%
%EXAMPLE
%   % no astigmatism
%   [phase amplitude]=tom_ctf2(-15e-6,0,0,5e-10,300000,[512 512],2e-3,.2,0.05,2.2e-3,0.8);
%   figure; tom_imagesc(phase);
%   figure; tom_imagesc(amplitude);
%
%   % with astigmatism
%   [phase amplitude]=tom_ctf2(-15e-6,5e-6,30,5e-10,300000,[512 512],2e-3,.2,0.05,2.2e-3,0.8);
%   figure; tom_imagesc(phase);
%   figure; tom_imagesc(amplitude);
%
%   % with astigmatism and decay - decay is not diectly applied!
%   [phase amplitude decay]=tom_ctf2(-15e-6,5e-6,30,5e-10,300000,[512 512],2e-3,.2,0.05,2.2e-3,0.8);
%   figure; tom_imagesc(phase);
%   figure; tom_imagesc(amplitude);
%   figure; tom_imagesc(decay);
%
%REFERENCES
% Reimer: Page 231, Frank: Three-dimensional electron microscopy of macromolecular assemblies,
% Edition: 2, illustrated, Oxford University Press US, 2006
% Page 36,37%
%
%SEE ALSO
%   TOM_CTF TOM_FOURIER TOM_IFOURIER TOM_CREATE_CTF
%
%
%   created by SN 04/06/09
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


voltage_relativistic=voltage*(1+voltage/1022000.0); % relativistic
lambda=sqrt(150.4/voltage_relativistic)*10^-10; % wavelength relativistic
Ny = 1./(2.0*pix_size);
Phi_0=Phi_0./180.*pi; % degree in rad

[x y] = meshgrid(-Ny:2*Ny/img_size(1):Ny-Ny/img_size(1), -Ny:2*Ny/img_size(2):Ny-Ny/img_size(2));
k = sqrt(x.^2.0+y.^2.0);
% Phase: Reimer, 6.41, page 218
% Amplitude: Reimer, 6.40
phase    =sin(pi/2.0.*(lambda.^3.*Cs.*k.^4 - 2.*lambda.*(Dz+(Dz_delta./2.0).*sin(2.0.*(atan2(x,y)-Phi_0))).*k.^2));
amplitude=cos(pi/2.0.*(lambda.^3.*Cs.*k.^4 - 2.*lambda.*(Dz+(Dz_delta./2.0).*sin(2.0.*(atan2(x,y)-Phi_0))).*k.^2));

% Frank: Three-dimensional electron microscopy of macromolecular
% assemblies, formula 2.4, 2.5 page 36. identical result !
%phase_f    =sin(pi.*2.0.*(lambda.^3.*Cs.*k.^4./4 - 0.5.*lambda.*(Dz+(Dz_delta./2.0).*sin(2.0.*(atan2(x,y)-Phi_0))).*k.^2));
%amplitude_f=cos(pi.*2.0.*(lambda.^3.*Cs.*k.^4./4 - 0.5.*lambda.*(Dz+(Dz_delta./2.0).*sin(2.0.*(atan2(x,y)-Phi_0))).*k.^2));
decay = 1;
if nargin>7
    decay = exp(-(k./(sigma*Ny)).^2);
end;
E_i_k=1;
if nargin>8
    % 'compound envelope function', Frank formula 2.14: E_k=E_i_k * E_e_k
    % q_0: a quantity of dimension 1/length specifying the size of the source
    % as it appears in the back focal plane.
    % partially coherent illumination, spatial coherence - assume 0.02mrad
    % q_0_r=0.01;
    q_0=q_0_r./((Cs.*lambda.^3).^0.25); % see formula 2.19a
    E_i_k=exp((-pi.^2.*q_0.^2.*(Cs.*lambda.^3.*k.^3-Dz.*lambda.*k).^2)); % for a gaussian source distribution
end;
% put in astigmatism !!!! here:
% E_i_k=exp((-pi.^2.*q_0.^2.*(Cs.*lambda.^3.*k.^3-Dz.*lambda.*k).^2)); % for a gaussian source distribution
% envelope due to energy spread, temporal coherence
%Cc=2.2e-3;
%deltaE=0.8;
E_k=1;
E_e_k=1;
if nargin>10
    delta_z = Cc.*deltaE./voltage;
    E_e_k=exp(-(pi.*delta_z.*lambda.*k.^2./2));
    E_k=E_i_k .* E_e_k;
end;
% all that in generalized coordinates - see Frank formula 2.18 ...
% Dz_h=Dz./((Cs.*lambda).^.5); % generalized defocus, see 2.18
% k_h=(Cs.*lambda.^3).^0.25.*k; % generalized spatial frequency, see 2.19
% q_0_r=0.01;
% q_0=q_0_r./((Cs.*lambda.^3).^0.25); % see formula 2.19a
% phase_h=sin(-pi.*Dz_h.*k_h.^2+pi./2.*k_h.^4);
% amplitude_h=cos(-pi.*Dz_h.*k_h.^2+pi./2.*k_h.^4);
% checked - is equal to phase!

%Amplitude in EM at an defocus of -15 mue                                                                        
%MEMSET 1000 READ Z (DIM 512) ERR 1 FOUR
%EMD (U 300 COE 2 VE 40000 L 10240) EMD (OBJ 5)
%SHIFT A B EMD B (DF -150000)
%PFKT B (TWOD) DIFF B C D (SPLIT)
%CROSS B COMPLETE B DIFF B C D (SPLIT)
%DIS (WIN 0 0 1)
%DIS C 0 0 1

%Amplitude in EM with an Dz_delta of 5 mue at an angle of 30 deg                                                                 
%MEMSET 1000 READ Z (DIM 512) ERR 1 FOUR
%EMD (U 300 COE 2 VE 40000 L 10240) EMD (OBJ 5)
%SHIFT A B EMD B (DF -150000 FA 50000 PHI 30)
%PFKT B (TWOD) DIFF B C D (SPLIT)
%CROSS B COMPLETE B DIFF B C D (SPLIT)
%DIS (WIN 0 0 1)
%DIS C 0 0 1

%the EM CTF and the TOM CTF deviates in the order of 10^-5. I assume that
%Matlab uses a higher numerical accuracy than EM or the sampling is slightly different. SN
% the term (Dz+(Dz_delta./2.0)... is defined as (Dz-(Dz_delta./2.0)... in
% EM - but Reimer and Frank have the TOM version.

% new function: include also envelope function:
% 1) spatial coherence
%  alpha: illumination aperture - assume 0.02mrad
%alpha=0.02e-3;
%Ks = exp(- ((pi*(Cs.*lambda.^2.*k.^3 - Dz.*k)*alpha).^2)/log(2));
% 2) temporal coherence
%deltaE=0.8; 
%Cc=2.2e-3;
%delta = Cc*deltaE/voltage;
%Kt = exp(- (pi*lambda*delta*k.^2/(4*log(2))).^2);
%K = Kt.*Ks;


