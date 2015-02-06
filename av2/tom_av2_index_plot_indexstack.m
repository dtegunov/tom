function tom_av2_index_plot_indexstack(index_stack,lookup,rows)
%TOM_AV2_INDEX_PLOT_INDEXSTACK plots index stack
%
%    tom_av2_index_plot_indexstack(index_stack,lookup,rows)
%PARAMETERS
%
%  INPUT
%   index_stack  index stack  by tom_av2_index_calc
%   lookup       lookup structure by tom_av2_index_calc
%   rows         vector with start and stop node [2 5]
%   
%  OUTPUT
%   -      
%
%EXAMPLE
%     tom_av2_index_plot_indexstack(stack,lookup,[2 5]);
%   
%
%REFERENCES
%
%SEE ALSO
%   tom_av2_index_search, tom_av2_index_calc 
%  
%
%   created by fb (eckster)
%   updated by ...
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

%determine size of matrix
base_max=max([lookup{:,3}]);
for i=1:size(lookup,1);
     tmp(i)=length(lookup{i,1});  
end;
max_x=max(tmp);
max_y=base_max^max_x;

x=rows(2)-rows(1)+1;
y=base_max^rows(2);


figure;
stack_count=0;
for i=1:rows(1)
    stack_count=stack_count+base_max^(i-1);
end;

row_count=rows(1);
row_count2=1;

for i=rows(1):rows(2)
    zz=((row_count2-1).*y)+round((y-base_max^row_count)./2)+1;    
    img_start=round((y-base_max^row_count)./2)+1;
    img_stop=round((y-base_max^row_count)./2)+base_max^i;
    for ii=img_start:img_stop
        subplot(x,y,zz); 
        try
            imagesc(index_stack(:,:,stack_count)'); 
        catch
            imagesc(zeros(size(index_stack,1),size(index_stack,2))); 
        end;
        
        axis image; colormap hot; set(gca,'xTick',[]); set(gca,'YTick',[]);
        numm=length(lookup{stack_count,2});
        title([lookup{stack_count,1} '||' num2str(numm)]);
        zz=zz+1;
        stack_count=stack_count+1;
    end
    row_count=row_count+1;
    row_count2=row_count2+1;
end;





