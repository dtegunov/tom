function tom_flicker_display(in1,in2,dsp_cub_flag,in_pause,num_iter,data_range,iso_view,iso_thr,iso_limvect,iso_camlight)
%tom_flicker_display displays images on the figure ...useful for small conf
%                     
%
%    tom_flicker_display(in1,in2,dsp_cub_flag,in_pause,num_iter)
%
%PARAMETERS
%
%  INPUT
%   in1            img or volume 1
%   in2            img or vlume 2
%   dsp_cub_flag   (0) 0,1,2 encodes the planes (...check tom_dspcub) 
%                   use -1 for isosurface representation 
%   in_pause       (0.2) pause after showing the image or volume in sec    
%   number_iter    (20) number of flicks  
%   data_range     (2std) for imagesc works only for 2d
%   iso_view       (-1) view for isosurface flicker -1 for interactive adjustment 
%   iso_thr        (mean + 2.2std) vector of thresh     
%   iso_limvect    (-1) limit axis (zoom) -1 for tight fit  
%   iso_camlight   ('headlight') 'right' , 'left' or light handle   
%
%
%  OUTPUT
%
%
%EXAMPLE
%
%%use dspcub 4 display
%tom_flicker_display(vol_1.Value,vol_2.Value,2,0.4,500)
%
%%use isosurface with interactive adj
%tom_flicker_display(vol_1.Value,vol_2.Value,-1,0.4,500,'',-1)
%
%use isosurface with fixed val
%tom_flicker_display(vol_1.Value,vol_2.Value,-1,0.4,'',500,[10 180],[2 2],[16.7026 49.4783 4.6319 61.0547 22.8909 46.6592],'right');
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by SN/FB 01/24/06
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

if (nargin < 3)
    dsp_cub_flag=0;
end;

if (nargin < 4)
   in_pause=0.2;
end;

if (nargin < 5)
    num_iter=20;
end;

if (nargin < 6)
    data_range='';
end;

if (nargin < 7)
    iso_view=-1;
end;

if (nargin < 8)
    iso_thr=((mean(in1(:))+2.5.*(std(in1(:)))) +  (mean(in2(:))+2.5.*(std(in2(:)))) )./2;
end;

if (nargin < 9)
    iso_limvect=[-1 -1 -1];
end;

if (nargin < 10)
    iso_camlight='headlight';
end;

if (length(iso_thr)==1)
    iso_thr=double([iso_thr iso_thr]);
end;

if (iso_view(1)==-1 || strcmp(iso_view,'spiral'))
    iso_view_vect=[0 0];
else
    iso_view_vect=iso_view;
end;


if (dsp_cub_flag==-1 )
    figure;
    [f1,v1]=isosurface(in1,iso_thr(1));
    [f2,v2]=isosurface(in2,iso_thr(2));
    ax_hand=iso_adj(f1,v1,in1,iso_view_vect,iso_camlight,[-1 -1 -1]);
    if (iso_view(1)==-1)
       R = input('Adjust Isosurface ...press return 2 continue','s');
    end;    
    x_lim=get(ax_hand,'Xlim');
    y_lim=get(ax_hand,'Ylim');
    z_lim=get(ax_hand,'Zlim');
    if (iso_limvect(1)==-1)
        iso_limvect=[x_lim y_lim z_lim];
    end;
    
    [v_a v_b]=view;
    iso_view_vect=[v_a v_b];  
    disp(['iso view: ' num2str([v_a v_b])]);
    disp(['iso thr: ' num2str(iso_thr)]);
    disp(['iso light: ' iso_camlight]);
    disp(['x y z lim: ' num2str(iso_limvect)]);
else
    figure;
    set(gcf,'position',[517   328   560   420]);
end;


if (iscell(in1)==0)
    
    if (min(size(in1))==1)
        dim_in='1d';
    else
        if (length(size(in1))==2)
            dim_in='2d';
        end;
        
        if (length(size(in1))==3)
            dim_in='3d';
        end;
    end;
 
    for i=1:num_iter
        
        if (strcmp(dim_in,'1d') && dsp_cub_flag > -1)
            plot(in1); drawnow; pause(in_pause); 
            plot(in2); drawnow; pause(in_pause); 
        end;
        
        if (strcmp(dim_in,'2d') && dsp_cub_flag > -1)
            if (isempty(data_range)) 
                tom_imagesc(in1); drawnow; pause(in_pause); 
                tom_imagesc(in2); drawnow; pause(in_pause); 
            else
                tom_imagesc(in1,'range',data_range); drawnow; pause(in_pause); 
                tom_imagesc(in2,'range',data_range); drawnow; pause(in_pause); 
            end;
        end;
        
        if (strcmp(dim_in,'3d') && dsp_cub_flag > -1)
            tom_dspcub(in1,dsp_cub_flag); drawnow; pause(in_pause); 
            tom_dspcub(in2,dsp_cub_flag); drawnow; pause(in_pause); 
        end;
       
        if (dsp_cub_flag ==-1)
            ax_hand=iso_adj(f1,v1,in1,iso_view_vect,iso_camlight,iso_limvect);   
            drawnow; pause(in_pause); delete(ax_hand);
            ax_hand=iso_adj(f2,v2,in1,iso_view_vect,iso_camlight,iso_limvect);
            drawnow; pause(in_pause); delete(ax_hand);
            if (strcmp(iso_view,'spiral'))
                if (mod(i,90)==0)
                    iso_view_vect=iso_view_vect+[40 0];
                end;
                if (mod(i,4)==0)
                    iso_view_vect=iso_view_vect+[0 20];
                end;
            end;
        end;
            
   end;
    
else
    %figure;
    if (min(size(in1{1}))==1)
        dim_in='1d';
    else
        if (length(size(in1{1}))==2)
            dim_in='2d';
        end;
        
        if (length(size(in1{1}))==3)
            dim_in='3d';
        end;
    end;
    
    for i=1:num_iter
        
        if (strcmp(dim_in,'2d') && dsp_cub_flag > -1)
            for ii=1:length(in1)
                if (isempty(data_range))
                    tom_imagesc(in1{ii}); drawnow; pause(in_pause);
                else
                    tom_imagesc(in1{ii},'range',data_range); drawnow; pause(in_pause);
                end;
            end;
        end;
        
        if (strcmp(dim_in,'3d') && dsp_cub_flag > -1)
            for ii=1:length(in1)
                tom_dspcub(in1{ii},dsp_cub_flag); drawnow; pause(in_pause);
            end;
        end;
    end;
    
end;


function ha=iso_adj(f,v,vol,iso_view_vect,iso_light,iso_lim)

p=patch('Faces',f,'Vertices',v);
ha=gca;
set(ha,'XTick',[]);
set(ha,'YTick',[]);
set(ha,'ZTick',[]);
set(p,'FaceColor','red','EdgeColor','none');
isonormals(vol,p);
daspect([1,1,1])
camlight('headlight');

view(iso_view_vect(1),iso_view_vect(2)); axis image;

lighting gouraud;
camlight(iso_light);
if (iso_lim(1) ~= -1)
    set(ha,'Xlim',iso_lim(1:2));
    set(ha,'YLim',iso_lim(3:4));
    set(ha,'ZLim',iso_lim(5:6));
end;



