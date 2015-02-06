function [particlepicker in_range]=tom_av2_ctf_fit_2_filefilter(dirpath,percentage_off)

%TOM_AV2_CTF_FIT_2_FILEFILTER creates a 'particlepicker' structure with
%.filelist and .filefilter based on the fitted defocus (by tom_fit_ctf_gui).
%
%   [particlepicker in_range]=tom_av2_ctf_fit_2_filefilter(dirpath,percentage_off)
%
%PARAMETERS
%
%  INPUT
%   dirpath             path of mircorgraphs, .em and .em.mat files
%   percentage_off      difference to 'inteded defocus' from the micrograph
%                       header in percent
%  
%  OUTPUT
%   particlepicker      structure with .filelist and .filefilter
%   in_range            total number of micrographs in range of threshold
%
%EXAMPLE
%   %create 'particlepicker' based on fitted defocus, with a difference less
%   %than 30 percent of the intended defocus value:
%
%   [particlepicker in_range]=tom_av2_ctf_fit_2_filefilter('070410/07042010/low/low_3*.em',30); 
%
%REFERENCES
%
%  ...
%
%SEE ALSO
%   tom_av2_filefilter_gui
%
%   created by SN/FB 04/13/10
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


d=dir(dirpath);

particlepicker.filelist={d.name};
[p n e]=fileparts(dirpath);
in_range=0;
for i=1:size(particlepicker.filelist,2)
    try
        load([p '/' particlepicker.filelist{i} '.mat']);
        
        Dz_fit=st_out.Fit.Dz_det;
        H=tom_reademheader([p '/' particlepicker.filelist{i}]);
        Dz_intended=H.Header.Defocus.*1e-10;
        if 100.*abs((Dz_intended-Dz_fit)./Dz_intended)<percentage_off
            particlepicker.filefilter{i}=[1];
            in_range=in_range+1;
        else
            particlepicker.filefilter{i}=[0];
        end;
    catch
        disp(['No tom_fit_ctf - defocus - .mat files. Already fitted?']);
    end;
end;


