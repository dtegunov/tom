function tom_makemovie3d(in,out,fps,quality,contrast,average,zoom,frames,direction,position,moviesize,bandpass,linewidth,redlinewidth,numbersflag)
%TOM_MAKEMOVIE3D generates an avi-file format movie file from a 3D array
%
%   tom_makemovie3d(in,out,fps,quality,contrast,average,zoom,frames,direction,position,moviesize,bandpass,linewidth,redlinewidth,numbersflag)
%
%   Scroll in z-direction and write out movie from single slices or 
%   average along the z-axis. Nice for visualizing 3D volumes of
%   a tomographic volume.
%
%PARAMETERS
%
%  INPUT
%   in                  ...
%   out                 ...
%   fps                 ...
%   quality             ...
%   contrast            ...
%   average             ...
%   zoom                ...
%   frames              ...
%   direction           ...
%   position            ...
%   moviesize           ...
%   bandpass            ...
%   linewidth           ...
%   redlinewidth        ...
%   numbersflag         ...
%  
%  OUTPUT
%    writes a movie directly to the harddrive
%
%EXAMPLE
%   tom_makemovie(in,'Movie.avi',10,100,[-1.5 1.5],3);
%
%   Make a movie 'Movie.avi' with 10 frames per second and 100 per cent
%   quality. Use -1.5 as a lower limit for the gray values and 1.5
%   as an upper limit and average 3 layers in z-direction
%
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by AK 05/13/05
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

outer_padding = 4;
inner_padding = 1;
outerdiff = 0;

%Prepare the figure
dimensions.x = size(in,1);
dimensions.y = size(in,2);
dimensions.z = size(in,3);

fig=figure;set(gcf,'Position',[50 50 moviesize(1) moviesize(2)]);
%set(fig,'DoubleBuffer','on');
set(fig,'Renderer','OpenGL'); 
%set(fig,'Renderer','zBuffer'); 
set(fig,'Color',[0 0 0]);
figuredimensions = get(gcf,'Position');

axesh.xy = axes();
axesh.xz = axes();
axesh.yz = axes();
axesh.box = axes();

figurewidth = figuredimensions(3) - outer_padding * 2 - inner_padding;
figureheight = figuredimensions(4) - outer_padding * 2 - inner_padding;
yfrac = dimensions.y / (dimensions.z + dimensions.y);
xfrac = dimensions.x / (dimensions.z + dimensions.x);
zfrac = dimensions.z / (dimensions.z + dimensions.y);

zsize = figureheight * zfrac;
xsize = figurewidth * xfrac;
ysize = figureheight * yfrac;


%Position of xy Slice
axes(axesh.xy);
axis ij;
left = outer_padding;
bottom = outer_padding+1 + zsize + inner_padding;
width = xsize;
height = ysize;
set(axesh.xy,'Units','pixel');
set(axesh.xy,'Tag','xy_axes');
set(axesh.xy, 'Position', [left bottom width height],'visible','on');
set(axesh.xy,'Units','normalized');
axis off;

%Position of xz slice
axes(axesh.xz);
axis ij;
left = outer_padding;
bottom = outer_padding+1;
width = xsize;
height = zsize;
set(axesh.xz,'Units','pixel');
set(axesh.xz, 'Tag', 'xz_axes');
set(axesh.xz, 'Position', [left bottom width height],'visible','on');
set(axesh.xz,'Units','normalized');
axis off;

%Position of yz slice
axes(axesh.yz);
axis ij;
left = outer_padding + inner_padding + xsize;
bottom = outer_padding+1 + zsize + inner_padding;
width = zsize;
height = ysize;
set(axesh.yz,'Units','pixel');
set(axesh.yz,'Tag','yz_axes');
set(axesh.yz, 'Position', [left bottom width height],'visible','on');
set(axesh.yz,'Units','normalized');
axis off;

%Position of 3D box
axes(axesh.box);
axis ij;
left = outer_padding + inner_padding + xsize;
bottom = outer_padding+1;
width = zsize;
height = zsize;
set(axesh.box,'Units','pixel');
set(axesh.box,'Tag','box');
set(axesh.box,'Position',[left bottom width height],'visible','on');
set(axesh.box,'Units','normalized');

%Prepare the avi handle
mov = avifile(out);
mov.Fps=fps;
if isequal(computer,'PCWIN')
mov.Compression='Cinepak';
else
mov.Compression='none';
end
% 'Cinepak', 'Indeo3', 'Indeo5', 'MSVC',', 'RLE', 'None'
mov.Quality=quality;


minframe = frames(1);
maxframe = frames(2);

if strcmp(direction,'x') == 1
    dir = 'yz';
    dir2 = 1;
    position.yz = minframe;
elseif strcmp(direction,'y') == 1
    dir = 'xz';
    dir2 = 2;
    position.xz = minframe;
else
    dir = 'xy';
    dir2 = 3;
    position.xy = minframe;
end

if minframe < 1 
    minframe = 1;
end

if maxframe > size(in,dir2) - average.(dir)
    maxframe = size(in,dir2) - average.(dir);
end

for lauf=minframe:maxframe
    render_slice(in,axesh,contrast,average.xz,'xz',position,zoom.xz,bandpass,linewidth,redlinewidth,numbersflag);
    render_slice(in,axesh,contrast,average.yz,'yz',position,zoom.yz,bandpass,linewidth,redlinewidth,numbersflag);
    render_slice(in,axesh,contrast,average.xy,'xy',position,zoom.xy,bandpass,linewidth,redlinewidth,numbersflag);
    F = getframe(gcf);
    mov = addframe(mov,F);
    drawnow;
    if strcmp(direction,'z')
        position.xy = position.xy+1;
    elseif strcmp(direction,'y')
        position.xz = position.xz+1;
    elseif strcmp(direction,'x')
        position.yz = position.yz+1;
    end
end

mov = close(mov);
close(fig);


%Render slices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function handles = render_slice(Volume,axesh,DataScale,sliceavg,orientation,position,actualaxis,bandpass,linewidth,redlinewidth,numbersflag)

%Volume:  image matrix
%axesh: axes handles
%DataScale: Contrast
%sliceavg: number of images to average 
%orientation: 'xy','xz','yz'
%position: struct xy, xz, yz
%bandpass: bandpass structure

persistent slice_xy;
persistent slice_yz;
persistent slice_xz;
persistent actualaxis_xy;
persistent actualaxis_yz;
persistent actualaxis_xz;

dimensions.x = double(size(Volume,1));
dimensions.y = double(size(Volume,2));
dimensions.z = double(size(Volume,3));


if orientation == 'xy'
		
	slice = squeeze(Volume(:,:,position.xy:position.xy+sliceavg));
	slice = mean(double(slice), 3);
	axes(axesh.xy);

	if bandpass.filter_xy ~= 0
		slice = tom_bandpass(double(slice), bandpass.filter_low, bandpass.filter_high);
		DataScale = DataScale ./ (dimensions.x*dimensions.y);
	end
	
	handles.images.xy = imagesc(slice',DataScale);colormap(gray);axis ij;
	set(axesh.xy,'xtick',[],'ytick',[]);
	axis(actualaxis);
	actualaxis_xy = actualaxis;
	slice_xy = tom_limit(slice',DataScale(1),DataScale(2));
	
	if numbersflag == 1
		h = findobj('Tag','text_xy2');
		delete(h);
		t = text(actualaxis(2)-(actualaxis(2)/100)*10,actualaxis(4)-(actualaxis(4)/100)*8,num2str(position.xy));
		set(t,'Tag','text_xy2','Color',[1 0 0],'FontWeight','bold','FontSize',16,'HorizontalAlignment','right');
	end
	
elseif orientation == 'xz'

	slice = Volume(:,position.xz:position.xz+sliceavg,:);
	slice = mean(double(slice), 2);
	slice = squeeze(slice);
	axes(axesh.xz);

	if bandpass.filter_xz ~= 0
			slice = tom_bandpass(double(slice), bandpass.filter_low, bandpass.filter_high);
			DataScale = DataScale ./ (dimensions.x*dimensions.z);
	end
	
	handles.images.xz = imagesc(slice',DataScale);colormap(gray);axis ij;
	set(axesh.xz,'xtick',[],'ytick',[]);	
	axis(actualaxis);
	actualaxis_xz = actualaxis;
	slice_xz = tom_limit(slice',DataScale(1),DataScale(2));

	if numbersflag == 1
		h = findobj('Tag','text_xz2');
		delete(h);
		t = text(actualaxis(2)-(actualaxis(2)/100)*10,actualaxis(4)-(actualaxis(4)/100)*10,num2str(position.xz));
		set(t,'Tag','text_xz2','Color',[1 0 0],'FontWeight','bold','FontSize',16,'HorizontalAlignment','right');
	end
	
elseif orientation == 'yz'

	
	slice = Volume(position.yz:position.yz+sliceavg,:,:);
	slice = mean(double(slice), 1);
	slice = squeeze(slice);
	axes(axesh.yz);
	
	if bandpass.filter_yz ~= 0
			slice = tom_bandpass(double(slice), bandpass.filter_low, bandpass.filter_high);
			DataScale = DataScale ./ (dimensions.y*dimensions.z);
	end
	
	handles.images.yz = imagesc(slice,DataScale);colormap(gray);axis ij;
	set(axesh.yz,'xtick',[],'ytick',[]);
	axis(actualaxis);
	actualaxis_yz = actualaxis;
	slice_yz = tom_limit(slice,DataScale(1),DataScale(2));

	if numbersflag == 1
		h = findobj('Tag','text_yz2');
		delete(h);
		t = text(actualaxis(2)-(actualaxis(2)/100)*12,actualaxis(4)-(actualaxis(4)/100)*8,num2str(position.yz));
		set(t,'Tag','text_yz2','Color',[1 0 0],'FontWeight','bold','FontSize',16,'HorizontalAlignment','right');
	end
	
end

%Draw Position markers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes(axesh.xy);
set(axesh.xy,'Units','pixel');
h = findobj('Tag','hline_xy2');
delete(h);
h = findobj('Tag','vline_xy2');
delete(h);

line('Xdata', [0 dimensions.x], 'Ydata', [position.xz position.xz], 'Color', 'r', 'LineWidth',redlinewidth, 'Tag', 'hline_xy2');
line('Xdata', [position.yz position.yz], 'Ydata', [0 dimensions.y], 'Color', 'r', 'LineWidth',redlinewidth, 'Tag', 'vline_xy2');

axes(axesh.xz);
set(axesh.xz,'Units','pixel');
h = findobj('Tag','hline_xz2');
delete(h);
h = findobj('Tag','vline_xz2');
delete(h);

line('Xdata', [0 dimensions.x], 'Ydata', [position.xy position.xy], 'Color', 'r', 'LineWidth',redlinewidth, 'Tag', 'hline_xz2');
line('Xdata', [position.yz position.yz], 'Ydata', [0 dimensions.z], 'Color', 'r', 'LineWidth', redlinewidth, 'Tag', 'vline_xz2');

axes(axesh.yz);
set(axesh.yz,'Units','pixel');
h = findobj('Tag','vline_yz2');
delete(h);
h = findobj('Tag','hline_yz2');
delete(h);

line('Xdata', [position.xy position.xy], 'Ydata', [0 dimensions.y], 'Color', 'r', 'LineWidth',redlinewidth, 'Tag', 'vline_yz2');
line('Xdata', [0 dimensions.z], 'Ydata', [position.xz position.xz], 'Color', 'r', 'LineWidth',redlinewidth, 'Tag', 'vline_yz2');

% Draw 3D Box
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if numel(actualaxis_xy) > 0 & numel(actualaxis_yz) > 0 & numel(actualaxis_xz) > 0 & orientation == 'xy'
	axes(axesh.box);
	newplot;
	hold on;
	
	%yz slice
	slice_yz = tom_mirror(slice_yz,'x');
	slice_yz = slice_yz';
	lauf = 0:5:dimensions.y;
	y = [lauf;lauf];
	x = [zeros(1,size(y,2))+position.yz;zeros(1,size(y,2))+position.yz];
	z = [zeros(1,size(y,2))+1;zeros(1,size(y,2))+dimensions.z];
	slice_yz(1:linewidth(1),:) = zeros()-1000;
	slice_yz(size(slice_yz,1)-linewidth(1):size(slice_yz,1),:)= zeros()-1000;
	slice_yz(:,size(slice_yz,2)-linewidth(1):size(slice_yz,2))= zeros()-1000;
	slice_yz(:,1:linewidth(1)) = zeros()-1000;
	warp(x,y,z,slice_yz);

	%xz slice
	lauf = 0:5:dimensions.x;
	x = [lauf;lauf];
	y = [zeros(1,size(x,2))+position.xz;zeros(1,size(x,2))+position.xz];
	z = [zeros(1,size(x,2))+1;zeros(1,size(x,2))+dimensions.z];

	pos_low = position.yz-linewidth(2);
	if pos_low < 1
		pos_low = 1;
	end
	pos_high = position.yz+linewidth(2);
	if pos_high > dimensions.x
		pos_high = dimensions.x
	end
	%slice_xz = tom_mirror(slice_xz,'y');
	slice_xz(:,pos_low:pos_high) = zeros()-1000;
	slice_xz(1:linewidth(1),:) = zeros()-1000;
	slice_xz(size(slice_xz,1)-linewidth(1):size(slice_xz,1),:)= zeros()-1000;
	slice_xz(:,size(slice_xz,2)-linewidth(1):size(slice_xz,2))= zeros()-1000;
	slice_xz(:,1:linewidth(1)) = zeros()-1000;
	warp(x,y,z,slice_xz);

	%xy slice
	lauf = 0:5:dimensions.x;
	x = [lauf;lauf];
	y = [zeros(1,size(x,2))+1;zeros(1,size(x,2))+dimensions.y];
	z = [zeros(1,size(x,2))+dimensions.z-position.xy;zeros(1,size(x,2))+dimensions.z-position.xy];
	pos_low = position.xz-linewidth(2);
	if pos_low < 1
		pos_low = 1;
	end
	pos_high = position.xz+linewidth(2);
	if pos_high > dimensions.y
		pos_high = dimensions.y;
	end
	pos_low2 = position.yz-linewidth(2);
	if pos_low2 < 1
		pos_low2 = 1;
	end
	pos_high2 = position.yz+linewidth(2);
	if pos_high > dimensions.x
		pos_high = dimensions.x;
	end
	%slice_xy = tom_mirror(slice_xy,'y');
	slice_xy(pos_low:pos_high,:) = zeros()-1000;
	slice_xy(:,pos_low2:pos_high2) = zeros()-1000;
	slice_xy(1:linewidth(1),:) = zeros()-1000;
	slice_xy(size(slice_xy,1)-linewidth(1):size(slice_xy,1),:)= zeros()-1000;
	slice_xy(:,size(slice_xy,2)-linewidth(1):size(slice_xy,2))= zeros()-1000;
	slice_xy(:,1:linewidth(1)) = zeros()-1000;
	warp(x,y,z,slice_xy);

	
	hold off;
	dx = actualaxis_xy(1);
	dy = actualaxis_xy(2);
	dz = actualaxis_xz(2);
	axis([1,dx+5,1,dy+5,1,dz+5]);axis off;
	set(gca,'Projection','orthographic');
	caxis(DataScale);
	%view(160,-30);
	axis tight;
	daspect([1,1,1]);
	%lighting phong;
	%camlight right;camlight left;camlight headlight;
end