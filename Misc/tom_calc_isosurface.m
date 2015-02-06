function [thresh,mass,vol_back,bin_mask]=tom_calc_isosurface(volume,molmass,pixelsize,accuracy,mask)
%TOM_CALC_ISOSURFACE search for a threshold including a specific molecular mass.
% Mass is 'white' !!!
%
%   [thresh,mass,vol_back,bin_mask]=tom_calc_isosurface(volume,molmass,pixelsize,accuracy,mask);
%
%PARAMETERS
%
%  INPUT
%   volume              3D Volume of particle
%   molmass             molecular mass in kDa
%   pixelsize           object pixelsize in Angstroem
%   accuracy            accuracy of search, good value is .1
%   mask                3D Mask Volume, threshold is calculated inside this
%                       mask volume only
%  
%  OUTPUT
%   thresh              calculated threshold
%   mass                calculated included mass
%   vol_back            thresholded volume
%   bin_mask            binarized, thresholded volume
%
%EXAMPLE
%   % load 20S proteasome x-ray structure:
%   % CRYSTAL STRUCTURE OF THE 20S PROTEASOME FROM YEAST AT 2.4 ANGSTROMS
%   % RESOLUTION, 15-Apr-1998, Groll, M., Ditzel, L., Loewe, J., Stock, D.,
%   % Bochtler, M., Bartunik, H.D., Huber, R.
%   proteasome = pdbread('http://www.rcsb.org/pdb/files/1ryp.pdb');
%   % proteasome = pdbread('http://www.rcsb.org/pdb/files/1pma.pdb');
%   % J. Loewe, Thermoplasma acidophilum 20S proteasome, 1995
%   % convert to electron density
%   emmap = tom_pdb2em2(proteasome, 2.1, 96); %Objectpixelsize at 2.1 A, Nyquist at 4.2 Angstrom
%   % filter to a resolution of 10A.
%   f=tom_norm(tom_bandpass(emmap,0,(96./2)./(10./4.2)),1);
%   [thresh mass vol vol_bin]=tom_calc_isosurface(f,720,2.1,.001);
%
%REFERENCES
%
%SEE ALSO
%   pdbread
%
%   created by SN 06/26/05
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

if nargin<5
    mask=ones(size(volume));
end;
volume_orig=volume;
volume=volume.*mask;
mass=max(max(max(volume)));
lauf=max(max(max(volume)));
while mass<molmass
    lauf=lauf-accuracy;
    nrmass=sum(sum(sum(volume>=lauf)));
    %    mass=nrmass;
    rho=1.3; % g/cm^3 for protein
    pixelsize_cm3=(pixelsize.*1.0e-10).^3./(0.01.^3);
    nr=nrmass;
    mass_g=nr.*rho.*pixelsize_cm3;
    mass_kg=mass_g./(10.^3);
    mass=mass_kg./(1.66e-27)./1000;
end;
thresh=lauf;
vol_back=(volume_orig>=thresh).*volume_orig;
bin_mask=(volume_orig>=thresh);
