function tom_av2_count_ppi(input_picklist,startingnumber,endingnumber)



%TOM_AV2_FIND_PAIRS finds pairs of particles based on a maximum distance in pixels.
%
%   tom_av2_count_ppi(input_picklist,startingnumber,endingnumber)
%
%PARAMETERS
%
%  INPUT
%   input_picklist                 input alignment file
%   output_picklist                output alignment file
%   max_dist                       maximum distance in pixels between two picked particles
%   new_raidus                     radius for new alignment file
%   show                           flag for displaying picked particles (1=yes, 0=no)
%
%EXAMPLE
%   tom_av2_count_ppi('../090902_p47f11/log/rec/15_corrf_high_128_clean_reorg.mat',3680,4566);
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
load(input_picklist);
counter=zeros(endingnumber,1);

for i=startingnumber:endingnumber;
    for k=zz+1:(zz+100);
        if k>size(align2d,2)
        break;
        end;
        
        if strcmp(align2d(1,k).filename,['/fs/sun17/lv01/pool/pool-nickell/26S/em/data/bohn/2d/100115_p47f11/high_corr/p47f11_' num2str(i) '.em'])==1;
            counter(i)=counter(i)+1;
            zz=k;
        end;
    end;
end;

for j=startingnumber:endingnumber;
    if counter(j)>70;
        disp([j counter(j)]);
    end;
end;

figure;plot(counter);
