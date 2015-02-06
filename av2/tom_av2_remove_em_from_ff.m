function tom_av2_remove_em_from_ff(flag,filefilter,proj_folder)
%TOM_AV2_REMOVE_EM_FROM_FF moves .em-files based on an input filefilter to a new folder.
%
%   tom_av2_remove_em_from_ff2(flag,filefilter,proj_folder)
%
%PARAMETERS
%
%  INPUT
%   flag                'singles' or 'pairs'
%   filefilter          path and name of input filefilter
%   proj_folder         project path  
% 
%EXAMPLE
%
%   tom_av2_remove_em_from_ff('singles','log/pick/08b_uncorr_low_128_509_ff.mat','/fs/pool/pool-nickell3/26S/em/data/eri/2d/100217_sc26s_rpn1gfp');
%
%REFERENCES
%
%SEE ALSO
%
%   created by SB 10/03/02
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

ff=load(filefilter);
 
target_high=[proj_folder '/high_del/'];
target_low=[proj_folder '/low_del/'];

disp([num2str(length(find([ff.particlepicker.filefilter{:}]==0))) ' bad images found!']);


for k=1:length(ff.particlepicker.filelist);
    if ff.particlepicker.filefilter{k}==0;
        if strcmp(flag,'pairs')
            source=[proj_folder '/high_corr/'  ff.particlepicker.filelist{k} '*'];
            call=['mv ' source ' ' target_high];
            disp(call);
            unix(call);
            source=[proj_folder '/high/'  ff.particlepicker.filelist{k} '*'];
            call=['mv ' source ' ' target_high];
            disp(call);
            unix(call);
        end;
        source=[proj_folder '/low_corr/'  ff.particlepicker.filelist{k} '*'];
        call=['mv ' source ' ' target_low];
        disp(call);
        unix(call);
        source=[proj_folder '/low/'  ff.particlepicker.filelist{k} '*'];
        call=['mv ' source ' ' target_low];
        disp(call);
        unix(call);
     end;
end;

disp([num2str(length(find([ff.particlepicker.filefilter{:}]==0))) ' bad images moved!']);