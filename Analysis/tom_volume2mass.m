function mass=tom_volume2mass(vol,pixelsize,thresh)
% tom_volume2mass calculates mass of given volume
%   
%      mass=tom_volume2mass(vol,pixelsize,thresh)
%
%  PARAMETERS
%  
%    INPUT
%     vol              input volume mass is white
%     pixelsize        pixs in Ang   
%     thresh           (0)  
%    
%    OUTPUT
%     mass            mass of the protein in kd
%
%  
%  EXAMPLE
%
%  mass=tom_volume2mass(tom_spheremask(ones(64,64,64),10),8.84);
%  
%  REFERENCES
%  
%  SEE ALSO
%     tom_calc_isosurface
%  
%     created by FB 04/11/2011
%  
%     Nickell et al., 'TOM software toolbox: acquisition and analysis for electron tomography',
%     Journal of Structural Biology, 149 (2005), 227-234.
%  
%     Copyright (c) 2004-2007
%     TOM toolbox for Electron Tomography
%     Max-Planck-Institute of Biochemistry
%     Dept. Molecular Structural Biology
%     82152 Martinsried, Germany
%     http://www.biochem.mpg.de/tom

if (nargin<3)
    thresh=0;
end;


nrmass=sum(sum(sum(vol>thresh)));
rho=1.3; % g/cm^3 for protein
pixelsize_cm3=(pixelsize.*1.0e-10).^3./(0.01.^3);
nr=nrmass;
mass_g=nr.*rho.*pixelsize_cm3;
mass_kg=mass_g./(10.^3);
mass=mass_kg./(1.66e-27)./1000;

