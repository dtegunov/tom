function [vol error_m]=tom_av2_backproj(align2d,count_st,disp_st)
%TOM_AV2_BACKPROJ creates ...
%
%   [vol error_m]=tom_av2_backproj(align2d,count_st,disp_st)
%
%PARAMETERS
%
%  INPUT
%   align2d             ...
%   count_st            ...
%   disp_st             ...
%  
%  OUTPUT
%   vol                 ...
%   error_m             ...
%
%EXAMPLE
%   ... = tom_av2_backproj(...);
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



%transfer often used data
error_m=0;
hist_count=count_st.hist;
h=tom_reademheader(align2d(hist_count,1).rec.file.Stack_Path);
size_stack=h.Header.Size;
h=tom_reademheader(align2d(hist_count,1).rec.file.Startmodel_Path);
size_vol=h.Header.Size;
angular_scan_eu=align2d(hist_count,1).rec.project.angular_scan_euler;
angular_scan=align2d(hist_count,1).rec.project.angular_scan;
weight_flag=align2d(hist_count,1).rec.backproj.weighting;
num_of_proj_per_class=align2d(hist_count,1).rec.classify.avg.num_array;
mask=tom_create_mask(align2d(hist_count,1).rec.backproj.filter.mask.mask);
min_num_of_proj=align2d(hist_count,1).rec.backproj.min_num_of_proj;
avg_path=align2d(hist_count,1).rec.classify.avg.path;
filter_st=align2d(hist_count,1).rec.backproj.filter.filter.filter;

vol=zeros(size_vol','single'); %allocate some memory 

%command window print
store.end=size(angular_scan,2); store.num_of_mesure=1;  store.mesured=0;
[store]=tom_disp_estimated_time(store,'start','backproj');

%display
disp_back_proj(1,'new',disp_st.developer);

%find used angleclasses ...build theta vector for weighting 
zz=1;
for i=1:size(angular_scan,2);
     if (num_of_proj_per_class(i) >= min_num_of_proj ) 
        eu_4_weight(zz,:)=angular_scan_eu(i,:);
        zz=zz+1;
     end;
end;



part=zeros(size_stack(1),size_stack(2));

% perform weighting and backprojection
for i=1:size(angular_scan,2)
    
        eu=align2d(hist_count,1).rec.project.angular_scan_euler(i,:);
        num_of_p=num_of_proj_per_class(i);

        if (num_of_p >= min_num_of_proj)    
        
            %read class avg
            part=tom_emreadc(avg_path,'subregion',[1 1 i],[size_stack(1)-1 size_stack(2)-1 0]);
            part=double(part.Value);
            part=(part./num_of_p); % weight it according to number
            if (strcmp(weight_flag,'Projection')==1)
                [w]=tom_calc_weight_functionc([size_stack(1) size_stack(2)],eu_4_weight,size_stack(1),eu);
                part=tom_apply_filter(part,filter_st);
                part=tom_apply_weight_function(part,w,ones(size(part)),0);
            end;
            part=single(part.*mask);
            vol=single(vol);
            tom_backproj3d_euler(vol,part,eu(1),eu(2),eu(3),[0 0 0]);

            %display
            disp_back_proj(part,'progress',disp_st.developer);
        end;

        
        %command window print
        store.i=i;
        store=tom_disp_estimated_time(store,'estimate_time');
        store=tom_disp_estimated_time(store,'progress');
    
end;

if (strcmp(weight_flag,'Volume')==1)
    [w]=tom_calc_weight_functionc(size_vol',eu_4_weight,size_stack(1));
    vol=tom_apply_filter(vol,filter_st);
    vol=tom_apply_weight_function(vol,w);
end;





function disp_back_proj(im,flag,disp_flag)

if (disp_flag==0)
    return;
end

if strcmp(flag,'new')
    if (isempty(findobj('tag','backproj')))
        figure; set(gcf,'tag','backproj'); 
        drawnow; set(gcf,'pos',[1235 143 565 434]); drawnow; drawnow;
    end;
 end;

if strcmp(flag,'progress')
    if (isempty(findobj('tag','backproj')))
        figure; set(gcf,'tag','backproj'); 
        drawnow; set(gcf,'pos',[1235 143 565 434]); drawnow; drawnow;
    end;
    figure((findobj('tag','backproj')) ); tom_imagesc(im,'noinfo');
end;

if strcmp(flag,'end')
    if (isempty(findobj('tag','backproj'))==0)
        %close(findobj('tag','backproj'));
    end;
end;

