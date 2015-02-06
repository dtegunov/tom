function cat_chimera_script(filename,models,translations,thresholds)

filename='./wow.py';
models{1}='/fs/sun11/lv01/pool/pool-nickell2/scratch/classify_from_blue_gene/ml3d_split_two_5/refine_bo/model_1/ProjMatch/run1/Iter_50/Iter_50_reconstruction.vol';
models{2}='/fs/sun11/lv01/pool/pool-nickell2/scratch/classify_from_blue_gene/ml3d_split_two_5/refine_bo/model_2/ProjMatch/run1/Iter_50/Iter_50_reconstruction.vol';
models{3}='/fs/sun11/lv01/pool/pool-nickell2/scratch/classify_from_blue_gene/ml3d_split_two_5/refine_bo/model_3/ProjMatch/run1/Iter_50/Iter_50_reconstruction.vol';
models{4}='/fs/sun11/lv01/pool/pool-nickell2/scratch/classify_from_blue_gene/ml3d_split_two_5/refine_bo/model_4/ProjMatch/run1/Iter_50/Iter_50_reconstruction.vol';
models{5}='/fs/sun11/lv01/pool/pool-nickell2/scratch/classify_from_blue_gene/ml3d_split_two_5/refine_bo/model_5/ProjMatch/run1/Iter_50/Iter_50_reconstruction.vol';

i=1;
for x=0:64:65.*5
translations.x(i)=x;
translations.y(i)=0;
translations.z(i)=0;
i=i+1;
end;


cmd=['cat /fs/pool/pool-bmsan-apps/tom_dev/Util/chimera_header.py > ' filename];
unix(cmd);

c=colormap('hsv');
idx_step=floor(size(c,1)./length(models));
idx=1;
for i=1:length(models)
    surface_colors.r(i)=c(idx,1);
    surface_colors.g(i)=c(idx,2);
    surface_colors.b(i)=c(idx,3);
    idx=idx+idx_step;
end;

for i=1:length(models)
    thresholds(i)=0.02;
end;

for i=1:length(models)
            call=['awk ' '''{' ...
              'gsub("File123.vol","' models{i} '");' ...  
              'gsub(" id_123 ","' num2str(i) '");' ...
              'gsub(" X_trans ","' num2str(translations.x(i)) '");' ...
              'gsub(" Y_trans ","' num2str(translations.y(i)) '");' ...
              'gsub(" Z_trans ","' num2str(translations.z(i)) '");' ...
              'gsub("surface_colors_r_replace","' num2str(surface_colors.r(i)) '");' ...
              'gsub("surface_colors_g_replace","' num2str(surface_colors.g(i)) '");' ...
              'gsub("surface_colors_b_replace","' num2str(surface_colors.b(i)) '");' ...
              'gsub("surface_levels_replace","' num2str(thresholds(i)) '");' ...
              'print }''' ' >> ' filename ' /fs/pool/pool-bmsan-apps/tom_dev/Util/chimera_volume.py'];

          unix(call);
end;
cmd=['cat /fs/pool/pool-bmsan-apps/tom_dev/Util/chimera_end.py >> ' filename];
unix(cmd);

unix('chimera --script wow.py &');
