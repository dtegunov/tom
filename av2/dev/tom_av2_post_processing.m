function [vol error_m]=tom_av2_post_processing(align2d,iter_num,vol);


if (exist('model','dir')==0)
    mkdir model;
end

symmetry=align2d(iter_num,1).model.sym;
sym_angle=align2d(iter_num,1).model.sym_angle;
mass=align2d(iter_num,1).model.mass;
pixel_size=align2d(iter_num,1).model.pixel_size;


mask_sp=tom_create_mask(align2d(iter_num,1).filter.mask.model_mask_sp);
mask_cy=tom_create_mask(align2d(iter_num,1).filter.mask.model_mask_cy);

file_path=align2d(iter_num,1).file_path;


error_m=0;

if (strcmp(symmetry,'C2')==1)
    vol=tom_symref(vol,symmetry,'C2',sym_angle);
end;


tom_emwrite([file_path '/model/model_' num2str(iter_num)],vol);

vol=vol.*mask_sp;
%[a b vol]=tom_calc_isosurface(-tom_norm(vol,1),mass,pixel_size,0.01);

[a b vol]=tom_calc_isosurface(tom_norm(-vol,1),mass,pixel_size,0.01);
%[a b vol]=tom_calc_isosurface(tom_norm(vol,1),mass,pixel_size,0.01);



%vol=vol.*tom_spheremask(ones(size(vol)),35,masksp_sigma);

vol=tom_rotate(mask_cy,[270 90 90]).*vol~=0.*mask_sp;
vol=-vol;
tom_norm(vol,1);
