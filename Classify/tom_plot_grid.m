function tom_plot_grid(pos,grid_size,grid_top,axes_handle, textflag, histogram)
%TOM_PLOT_GRID plots a 2d grid
%
%   tom_plot_grid(pos,axes_handle,grid_size,grid_top)
%
%   plots positions as a 2d grid
%   useful for ploting a som
%
%PARAMETERS
%
%  INPUT
%   radius_in     a scalar defining the inner radius
%   radius_out    a scalar defining the outer radius
%   overall_size  size of 3D array
%  
%  OUTPUT
%   cyl           a 3D array with the cylinder
%
%EXAMPLE
%   cyl=tom_cylinder(8, 10, [32 32 32]);
%   creates a 32 cube cyl with a cylinder of an inner
%   radius of 8 and an outer radius of 10.
%
%REFERENCES
%
%SEE ALSO
%   TOM_CIRCLE
%
%   created by SN 04/04/04
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

%generate colormap
figure;cm=colormap('hot');close(gcf);

if nargin < 5
    textflag = 1;
end

if (nargin < 4)
    axes_handle=gca;
end;

if (nargin < 3)
    grid_top='rect';
end;

if (nargin < 2)
    grid_size=sqrt(size(all_pos,1));
end;

if nargin < 6
    histogram = ones(grid_size(1).*grid_size(2),1);
end

histogram = round(tom_norm(histogram,63)+1);
axes(axes_handle);

hold on;
for i=1:size(pos,1)
    tmp_neig=tom_neighbour2dgrid(i,grid_size,grid_top);
    x_tmp(1)=pos(i,1);
    y_tmp(1)=pos(i,2);
    for ii=1:max(size(tmp_neig))
        x_tmp(2)=pos(tmp_neig(ii),1);
        y_tmp(2)=pos(tmp_neig(ii),2);
      plot(x_tmp,y_tmp,'LineWidth',1.5);
    end;
end

for i=1:size(pos,1)
    tmp_neig=tom_neighbour2dgrid(i,grid_size,grid_top);
    x_tmp(1)=pos(i,1);
    y_tmp(1)=pos(i,2);
    for ii=1:max(size(tmp_neig))
        x_tmp(2)=pos(tmp_neig(ii),1);
        y_tmp(2)=pos(tmp_neig(ii),2);
     % plot(x_tmp,y_tmp,'LineWidth',1.5);
    end;

    plot(pos(i,1),pos(i,2),'ro','MarkerSize',8,'MarkerFaceColor',cm(histogram(i),:),'MarkerEdgeColor','b');
    if textflag == 1
        text(pos(i,1)+abs(max(pos(:,1)).*0.05),pos(i,2)-abs(max(pos(:,2)).*0.05),['(' num2str(floor((i-1)./grid_size(1))+1) ',' num2str( i - floor((i-1)./grid_size(1)).*grid_size(1) ) ')']);
    end
    
 end;
hold off;
