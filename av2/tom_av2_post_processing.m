function [vol error_m]=tom_av2_post_processing(align2d,count_st,vol);
%TOM_AV2_POST_PROCESSING creates ...
%
%   [vol error_m]=tom_av2_post_processing(align2d,count_st,vol);
%
%PARAMETERS
%
%  INPUT
%   align2d             ...
%   count_st            ...
%   vol                 ...
%  
%  OUTPUT
%   vol                 ...
%   error_m             ...
%
%EXAMPLE
%   ... = tom_av2_post_processing(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by .. mm/dd/yy
%   updated by ..
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

%transfer variables
hist_count=count_st.hist;
output_dir=align2d(hist_count,1).rec.file.Outputdir;
sym_flag=align2d(hist_count,1).rec.post_proc.sym.on;
sym=align2d(hist_count,1).rec.post_proc.sym.symmetry;
sym_angle=align2d(hist_count,1).rec.post_proc.sym.rot_angle;
bin_flag=align2d(hist_count,1).rec.post_proc.binarize.on;
bin_mass=align2d(hist_count,1).rec.post_proc.binarize.mass;
bin_pixel_size=align2d(hist_count,1).rec.post_proc.binarize.pixelsize;
bin_mass_is_black=align2d(hist_count,1).rec.post_proc.binarize.mass_is_black;
sign_f=1;
error_m=0;


mask=tom_create_mask(align2d(hist_count,1).rec.post_proc.filter.mask);

tom_emwrite([output_dir '/step' num2str(count_st.step) '/model/model_' num2str(count_st.iteration) '.em' ],vol);



% sym model
if (sym_flag)
    vol=tom_rotate(vol,sym_angle);
    vol=tom_symref(vol,sym);
    vol=tom_rotate(vol,[-sym_angle(2) -sym_angle(1) -sym_angle(3)]);
end;

vol=vol.*mask;


% binarize model
if (bin_flag)
    if (bin_mass_is_black==1)
        sign_f=-1;
    end;
    [a b vol]=tom_calc_isosurface(tom_norm((vol.*sign_f),1),bin_mass,bin_pixel_size,0.01);
    vol=(vol~=0);
    if (bin_mass_is_black==1)
        vol=-vol;
    end;
    tom_norm(vol,1);
end;

tom_emwrite([output_dir '/step' num2str(count_st.step) '/model/model_bin_' num2str(count_st.iteration) '.em' ],vol);



[a,b]=system(['chmod -R ugo+rwx ' output_dir '/step' num2str(count_st.step) '/model/']);