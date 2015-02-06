function tom_av2_clean_picklist(input_picklist,output_picklist)

%TOM_AV2_CLEAN_PICKLIST removes particles from alignment file, if
%image doesn't exist.
%
%   tom_av2_clean_picklist(input_picklist,output_picklist)
%
%PARAMETERS
%
%  INPUT
%   input_picklist        input alignment file
%   output_picklist       output alignment file
%
%
%EXAMPLE
%   tom_av2_clean_picklist('25_corrf_high_128.mat','25_corrf_high_128_clean.mat');
%
%REFERENCES
%
%SEE ALSO
%
%   created by SB 10/02/24
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

zz=0;
input_picklist=load(input_picklist);
align2e=input_picklist.align2d;

for k=1:size(input_picklist.align2d,2);
    [a b]=unix(['ls ' input_picklist.align2d(1,k).filename]);
    if a==0;
        zz=zz+1;
        align2e(1,zz)=input_picklist.align2d(1,k);
  
    else
      disp(['Image ' input_picklist.align2d(1,k).filename ' not existent.']); 
    end;
    if mod(k,1000)==0;
        disp(['Particle ' num2str(k) ]);
    end;
end;

disp(['Particles deleted:' num2str(k-zz)]);
tmp=align2e(1,1:zz);
align2e=tmp;
align2d=align2e;
save(output_picklist,'align2d');