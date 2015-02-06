function final_mask=tom_mass_grow(vol,seed_mask,final_mass,thr,pixs)
%tom_mass_grow 
%   [vol_clean dust]=tom_mass_grow(vol,mass,bw_cer,max_iter,verbose)
%PARAMETERS
%
%  INPUT
%   seed               input volume or image
%   final_mass         mass in kd
%   thr                threshold
%   pixs               pixelsize in Ang
%
% 
%  OUTPUT 
%   final_mask         mask
%
%
%
%EXAMPLE
%   
% final_mask=tom_mass_grow(vol_1.Value,vol_clean,270,0.016,2.21);
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by FB late june 10
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

rho=1.3; % g/cm^3 for protein

ini_grow_th=0.13;
max_mass_grow=10;

min_mass_grow=2;

grow_th=ini_grow_th;
mass_old=0;
mask_old=seed_mask;

for i=1:5
    
    fprintf('      %s','Adapting mass grow...');
    for ii=1:60
        mask=tom_filter(mask_old,4);
        mask=mask>grow_th;
        vol_bin=(vol.*mask)>thr;
        mask=vol_bin;
        num_vox=length(find(vol_bin));
        
        pixelsize_cm3=(pixs.*1.0e-10).^3./(0.01.^3);
        mass_g=num_vox.*rho.*pixelsize_cm3;
        mass_kg=mass_g./(10.^3);
        mass=mass_kg./(1.66e-27)./1000;
        
        fprintf(' %3.1f  ',(mass-mass_old));
        if ( ((mass-mass_old) > min_mass_grow) &&  ((mass-mass_old) < max_mass_grow ))
            fprintf(' done \n');
            break;
        end;
        
        if ((mass-mass_old)>max_mass_grow)
           % disp(['  mass grow: ' num2str(mass-mass_old) ]);
            grow_th=(grow_th.*1.07);
            
        end;
        if ((mass-mass_old) < min_mass_grow)
            %disp(['  mass grow: ' num2str(mass-mass_old) ]);
            grow_th=(grow_th.*0.93);
          
        end;
    end;
    if (ii==40)
        warning('Mass grow could not be adapted');
    end;
    
    mass_old=mass;
    mask_old=mask;
    
    disp(['Mass: ' num2str(mass)  ' num of voxels: ' num2str(num_vox) ]);
    if (mass > final_mass)
        break;
    end;
    
end;
final_mask=vol_bin;


