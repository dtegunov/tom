function corrected=tom_correct_for_ctf_and_mtf(in,defocus,objectpixelsize,microscope,cutoff)

%TOM_CORREECT_FOR_CTF_AND_MTF  Corrects a volume for ctf and mtf
%
%   corrected=tom_correct_for_ctf_and_mtf(in,defocus,microscope,cutoff)
%
%
%PARAMETERS
%
%  INPUT
%   in                  original volume, only cubic
%   defocus             defocus of the reconstruction in mue, underfocus is
%                       negative
%   objectpixelsize     Objectpixelsize of reconstruction in nm
%   microscpe           electron microscope type, only 'Polara' implemented
%                       yet
%   cutoff              cutoff frequency in pixel
%
%  OUTPUT
%   corrected           corrected volume
%
%EXAMPLE
%       corrected=tom_correct_for_ctf_and_mtf(in,-4,0.72,'Polara',20)
%
%
%REFERENCES
%
%SEE ALSO
%
%   created by SN 16/01/07
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

sx=size(in,1);

if isequal(microscope,'Polara')
%    mtf_fct=tom_emread('polara_mtf.em');
    %load polara_mtf
    load mtf_Eagle_SN109_200kV
    mtf=mtf_Eagle_SN109_200kV;
    %mtf = mtf';
    idx = find(isnan(mtf));
    mtf(idx) = 0;
    mtf_fct=imresize(mtf,((sx./2)./size(mtf,1)));
    mtf = zeros((sx./2),2.*sx,sx);


    for ii=1:2.*sx 
        for jj=1:sx
            mtf(:,ii,jj) = mtf_fct; 
        end;
    end;




    mtf_cart = tom_sph2cart(mtf);


    %ctf_fct = tom_ctf(defocus,objectpixelsize, 300,sx,2);
%    ctf_fct = tom_ctf(defocus,objectpixelsize, 160,sx,2);
    %ctf_fct=imresize(ctf_theo,(40./1024));
    %ctf = zeros(sx./2,2.*sx,sx);

    %ctf_fct=tom_filter(abs(ctf_fct),3,'quadr','real');


    %q=(ctf_fct(1,1:sx./2)./(abs(ctf_fct(1,1:sx./2))+10e-6 ));
    %q(1)=1;

    %for ii=1:2.*80, for jj=1:80, ctf(:,ii,jj) = abs(ctf_fct(1,1:40)); end;end;
    %for ii=1:2.*sx, for jj=1:sx, ctf(:,ii,jj) = q; end;end;




  %  ctf_cart = tom_sph2cart(ctf);

    disp('mtf only !!!');
    all=1./(mtf_cart+0.000001);

%    all=(1./(mtf_cart+0.000001));

    mask=tom_spheremask(ones(sx,sx,sx),cutoff,1);
    all_band=all.*mask;

%    in=tom_emread('26S.vol');
    corrected=tom_apply_weight_function(in,all_band);



else
    error(['Electron microscope not known.']);
end;

