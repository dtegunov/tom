function tom_ctf_autofit_correctphase_swedish(directory,min_freq,max_freq,enhance_filter_min,enhance_filter_max,enhance_weight,outdir,method,epsilon)
% TOM_CTF_AUTOFIT searches for em images in a given directory and tries to
% do a CTF fit on each image. If successful the new defocus value is
% written to the header field "FocusIncrement", on failure the nominal
% value is copied.
%
% tom_ctf_autofit(directory, threshold, filterval)
%
%
%  INPUT
%
%  directory    The full path to the directory containing the images
%  threshold    The maximum difference between fitted and nominal defocus
%               in mu (optional) (default: 2)
%  filterval    The kernel size for a real space quadratic kernel to be
%               applied to each image before CTF fitting (optional) (default 1)
%  startval     start defocus value in nm (optional)
%
%  OUTPUT
%
%EXAMPLE
%   tom_ctf_autofit(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by AK 09/11/06
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




dircell = tom_HT_getdircontents(directory,{'em'});


if nargin < 9
    epsilon = 0;
end
if nargin < 8
    method = 'leave';
end

parfor (i=1:length(dircell))
    
    %if tom_isemfile(['/fs/sun17/lv01/pool/pool-nickell/177-36mer/180208-210208/low_corr/test_' num2str(i) '.em']) ~= 1
    %    continue;
    %end
    
    if i == 822
        continue;
    end
    
    file = [directory '/test_' num2str(i) '.em'];
    %file = dircell{i};
    if exist(['/fs/sun17/lv01/pool/pool-nickell/177-36mer/180208-210208/low_corr/test_' num2str(i) '.em'],'file')
    	continue;
    end

    disp(file);
    header = tom_emreadc(file);
    %header = tom_emreadc([directory '/' file]);
    Dz = header.Header.Defocus;
    voltage = header.Header.Voltage./1000;
    objectpixelsize = header.Header.Objectpixelsize;
    Cs = header.Header.Cs;
    Ca = 2;
    ctfmodelsize = 0;

    header.Value = tom_xraycorrect(header.Value);
    
    st = tom_xmipp_adjust_ctf(tom_calc_periodogram(header.Value,256),Dz,voltage,objectpixelsize,ctfmodelsize,Cs,min_freq,max_freq,Ca,enhance_filter_min,enhance_filter_max,enhance_weight);
    
%     if i<692
%         st.DeltafU = st.DeltafU + 20000;
%         st.DeltafV = st.DeltafV + 20000;
%     elseif i<1070
%         st.DeltafU = st.DeltafU + 30000;
%         st.DeltafV = st.DeltafV + 30000;
%     else
%         st.DeltafU = st.DeltafU + 50000;
%         st.DeltafV = st.DeltafV + 50000;
%     end
    if i<440
        st.DeltafU = st.DeltafU + 20000;
        st.DeltafV = st.DeltafV + 20000;
    elseif i<791
        st.DeltafU = st.DeltafU + 30000;
        st.DeltafV = st.DeltafV + 30000;
    elseif i<1175
        st.DeltafU = st.DeltafU + 20000;
        st.DeltafV = st.DeltafV + 20000;
    else
        st.DeltafU = st.DeltafU + 50000;
        st.DeltafV = st.DeltafV + 50000;
    end
        


    header = tom_emreadc(['/fs/sun17/lv01/pool/pool-nickell/177-36mer/180208-210208/low/test_' num2str(i) '.em']);
    
    im_corr = tom_xmipp_ctf_correct_phase(header.Value,st,method,epsilon);
    header.Value = tom_xraycorrect(header.Value);
    tom_emwrite(['/fs/sun17/lv01/pool/pool-nickell/177-36mer/180208-210208/low_corr/test_' num2str(i) '.em'],im_corr);
    %tom_emwrite([outdir '/' file],im_corr);
    disp([outdir '/' file]);
    
end





