function tom_chimera_volumes(filename,models,translations,thresholds,colormap_models,bright_trans,origins,hide_dust,bild_file)
%TOM_CHIMERA_VOLUMES creates a Chimera .py-script for displaying volumes and calls
%  Chimera
%
%   tom_chimera_volumes(filename,models,translations,thresholds,colormap_models,bright_trans,origins,hide_dust,bild_file)
%
%
%  INPUT
%   filename:            output file name of Chimera .py-script
%   models:              cell with filenames of 3D volumes. Volume data in
%                        EM or Spider format!!!
%   translations:        translations as 3d vector or 'linear_x', 'linear_y', 'linear_z'.
%   thresholds:          thresholds values as vector
%   colormap_models:     colormap for models, eg. 'hsv', 'hot' ...
%   bright_trans         (opt.) brightness and transperency level as vector [0.5 0.3] 
%   origins              (opt.) origions of volums   
%   hide_dust            (opt.) hide_dust values as vector
%   bild_file            (opt.) build file to overlay with structures      
%
%  OUTPUT
%   creates the Chimera .py-script and calls Chimera
%
%EXAMPLE
%
% models{1}='/fs/sun05/lv01/pool/bmsan/apps/tom_dev/Chimera/test_files/model_1.spi';
% models{2}='/fs/sun05/lv01/pool/bmsan/apps/tom_dev/Chimera/test_files/model_2.spi';
% models{3}='/fs/sun05/lv01/pool/bmsan/apps/tom_dev/Chimera/test_files/model_3.spi';
% models{4}='/fs/sun05/lv01/pool/bmsan/apps/tom_dev/Chimera/test_files/model_1.em';
%
% translations='linear_x'
%
% thresholds=[0.022]
% colormap_models='hsv'
%
% tom_chimera_volumes('show.py',models,translations,thresholds,colormap_models)
%
% % or, the same volume at different thresholds:
%
% models{1}='/fs/sun05/lv01/pool/bmsan/apps/tom_dev/Chimera/test_files/model_1.em';
% models{2}='/fs/sun05/lv01/pool/bmsan/apps/tom_dev/Chimera/test_files/model_1.em';
% models{3}='/fs/sun05/lv01/pool/bmsan/apps/tom_dev/Chimera/test_files/model_1.em';
% models{4}='/fs/sun05/lv01/pool/bmsan/apps/tom_dev/Chimera/test_files/model_1.em';
%
% translations='linear_y'
%
% thresholds=[0.02 0.04 0.08 0.1]
% colormap_models='jet'
%
% tom_chimera_volumes('show.py',models,translations,thresholds,colormap_models)
%
% Gold:
% colormap_models=[0.8 0.64 0.2]
%
%REFERENCES
%
%SEE ALSO
%     
%
%   created by SN,FF 12/02/10
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

if (isempty(filename))
    filename='show.py';
end;

if tom_isspiderfile(models{1})
    in=tom_spiderread(models{1});
end;

if tom_isemfile(models{1})
    in=tom_emread(models{1});
end;

if tom_iseman2h5file(models{1})
   in=tom_eman2_read(models{1}); 
end;

if (nargin <3)
    translations='';
end;

if (nargin <4)
    thresholds=2.5;
end;


if (nargin<5)
    colormap_models=[1 0 0];
end;

if (nargin<6)
    bright_trans='';
end;

if (nargin<7)
    origins='';
end;

if (nargin<8)
    hide_dust='';
end;

if (nargin<9)
    bild_file='';
end;



if (isempty(bright_trans))
    clear('bright_trans')
    for i=1:length(models)
        bright_trans(1,i)=1;
        bright_trans(2,i)=0;
    end; 
end;

if (isempty(origins))
    clear('origins')
    for it=1:length(models)
        origins(1,it)=0;
        origins(2,it)=0;
        origins(3,it)=0;
     end;
end;

if (isempty(translations))
    clear('translations');
    for itrt=1:length(models)
        translations(1,itrt)=0;
        translations(2,itrt)=0;
        translations(3,itrt)=0;
     end;
end;

if isnumeric(translations(1,1))==0
    if strcmp(translations,'linear_x')
        translations=0;
        i=1;
        for x=0:in.Header.Size(1)./2:(in.Header.Size(1)./2+1).*size(models,2)
            translations(1,i)=x;
            translations(2,i)=0;
            translations(3,i)=0;
            i=i+1;
        end;
    end;
    if strcmp(translations,'linear_y')
        translations=0;
        i=1;
        for x=0:in.Header.Size(1)./2:(in.Header.Size(1)./2+1).*size(models,2)
            translations(1,i)=0;
            translations(2,i)=x;
            translations(3,i)=0;
            i=i+1;
        end;
    end;
    if strcmp(translations,'linear_z')
        translations=0;
        i=1;
        for x=0:in.Header.Size(1)./2:(in.Header.Size(1)./2+1).*size(models,2)
            translations(1,i)=0;
            translations(2,i)=0;
            translations(3,i)=x;
            i=i+1;
        end;
    end;
end;



head_py_scrt='/fs/pool/pool-bmsan-apps/tom_dev/Chimera/py_files/chimera_header.py';
vol_head_py_scrt='/fs/pool/pool-bmsan-apps/tom_dev/Chimera/py_files/chimera_vol_header.py';
vol_py_scrt='/fs/pool/pool-bmsan-apps/tom_dev/Chimera/py_files/chimera_volume.py';
end_py_scrt='/fs/pool/pool-bmsan-apps/tom_dev/Chimera/py_files/chimera_end.py';


cmd=['cat ' head_py_scrt ' > ' filename];
unix(cmd);

%add restore 4 hide dust
if (isempty(hide_dust)==0)
    add_hide_dust_mod(filename,hide_dust)
end;


if size(thresholds,2)~=size(models,2)
    disp('assume first given threshold as threshold for all volumes');
    thresholds=repmat(thresholds(1),[size(models,2) 1]);
end;

c=colormap(colormap_models);
try
    close(gcf);
catch Me
end;

idx_step=floor(size(c,1)./length(models));
idx=1;
for i=1:length(models)
    surface_colors.r(i)=c(idx,1);
    surface_colors.g(i)=c(idx,2);
    surface_colors.b(i)=c(idx,3);
    idx=idx+idx_step;
end;

cmd=['cat ' vol_head_py_scrt ' >> ' filename];
unix(cmd);


for i=1:length(models)
    if tom_isspiderfile(models{i})
        in=tom_spiderread(models{i});
    end;
    if tom_isemfile(models{i})
        in=tom_emread(models{i});
    end;
    if tom_iseman2h5file(models{i})
        in=tom_eman2_read(models{i});
    end;
    [a b c]=fileparts(models{i});
    
    call=['awk ' '''{' ...
        'gsub("File123.vol","' models{i} '");' ...
        'gsub(" id_123 ","' num2str(i) '");' ...
        'gsub("->Name123.vol<-","' [b c] '");' ...
        'gsub(" X_trans ","' num2str(translations(1,i)) '");' ...
        'gsub(" Y_trans ","' num2str(translations(2,i)) '");' ...
        'gsub(" Z_trans ","' num2str(translations(3,i)) '");' ...
        'gsub(" X_origin ","' num2str(origins(1,i)) '");' ...
        'gsub(" Y_origin ","' num2str(origins(2,i)) '");' ...
        'gsub(" Z_origin ","' num2str(origins(3,i)) '");' ...
        'gsub("surface_colors_r_replace","' num2str(surface_colors.r(i)) '");' ...
        'gsub("surface_colors_g_replace","' num2str(surface_colors.g(i)) '");' ...
        'gsub("surface_colors_b_replace","' num2str(surface_colors.b(i)) '");' ...
        'gsub("surface_levels_replace","' num2str(thresholds(i)) '");' ...
        'gsub("surface_brightness_factor_replace","' num2str(bright_trans(1,i)) '");' ...
        'gsub("transparency_factor_replace","' num2str(bright_trans(2,i)) '");' ...
        'gsub("SZ_m1_mod","' num2str(in.Header.Size(1)-1) '");' ];
    
    if tom_isemfile(models{i})
        call=[call 'gsub("spider","tom_em");'];
        call=[call 'gsub("my_grid_id_123","");'];
    end;
    
    if tom_isspiderfile(models{i})
        call=[call 'gsub("spider","spider");'];
        call=[call 'gsub("my_grid_id_123","");'];
    end;
    
    if tom_iseman2h5file(models{i})
        call=[call 'gsub("spider","emanhdf");'];
        call=[call 'gsub("my_grid_id_123","/MDF/images/0/image");'];
    end;
    
    call=[call 'print }''' ' >> ' filename ' ' vol_py_scrt];
    unix(call);
    
end;
cmd=['cat ' end_py_scrt ' >> ' filename];
unix(cmd);

cmd='chimera ';
if (isempty(bild_file))
    cmd=[cmd ' --script ' filename ' &'];
else
    if (iscell(bild_file))
        cmd=[cmd filename];
        for i=1:length(bild_file)
            cmd=[cmd ' ' bild_file{i}];
        end;
        cmd=[cmd ' &'];
    else
        cmd=[cmd filename ' ' bild_file ' &'];
    end;
end;

disp(cmd);
unix(cmd);

function add_hide_dust_mod(filename,hide_dust)

fid=fopen(filename,'a');
fprintf(fid,'%s\n','def restore_hide_dust():');
fprintf(fid,' %s\n','hide_dust_state = \');
fprintf(fid,'    %s\n','{');
fprintf(fid,'     %s,\n','''class'': ''Hide_Dust_State''');
fprintf(fid,'     %s\n','''dust_table'': {');
for i=1:length(hide_dust)
    fprintf(fid,['        ( %d, 0, ): ( ' '' '''size''' '' ', %f, ),' ],i,hide_dust(i));
end;
fprintf(fid,'     %s\n','},');
fprintf(fid,'     %s\n',['''version''' ': 1,']);
fprintf(fid,'    %s\n','}');
fprintf(fid,' %s\n','try:');
fprintf(fid,'  %s\n','import HideDust.session');
fprintf(fid,'  %s\n','HideDust.session.restore_hide_dust_state(hide_dust_state)');
fprintf(fid,' %s\n','except:');
fprintf(fid,'  %s\n',['reportRestoreError(''' 'Error restoring hide dust''' ')']);
fprintf(fid,'%s\n',[]);
fprintf(fid,'%s\n','registerAfterModelsCB(restore_hide_dust)');
fclose(fid);

