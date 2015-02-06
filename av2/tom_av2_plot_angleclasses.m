function [st_tmp cccs]=tom_av2_plot_agleclasses(align2d,flag,iter_num,class_nr,num_of_parts)
%TOM_AV2_PLOT_ANGLECLSSES creates ...
%
%   [st_tmp cccs]=tom_av2_plot_agleclasses(align2d,flag,iter_num,class_nr,num_of_parts)
%
%PARAMETERS
%
%  INPUT
%   align2d             ...
%   flag                ...
%   iter_num            ...
%   class_nr            ...
%   num_of_parts        ...
%  
%  OUTPUT
%   st_tmp              ...
%   cccs                ...
%
%EXAMPLE
%   ... = tom_av2_plot_agleclasses(...);
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

st_tmp=0;
cccs=0;

if (nargin==2)
    flag='';
end;


if (strcmp(flag,'show_all_classes'))
    figure;
    for i=iter_num(1):iter_num(2)
        im=tom_emreadc(align2d(i,1).rec.classify.avg.path);
        for iii=1:size(im.Value,3)
            num_of_p=align2d(i,1).rec.classify.avg.num_array(iii);
             im.Value(:,:,iii)=im.Value(:,:,iii)./num_of_p;
        end;
        tom_dspcub(im.Value);
        drawnow;
        pause(0.8);
        set(gcf,'name',num2str(i))
    end;
end;



if (strcmp(flag,'show_all_classes_dist'))
    for i=iter_num(1):iter_num(2)
        num=align2d(i,1).rec.classify.avg.num_array;
        hold on;
        plot(num);
        hold off;
    end;
end;



if (strcmp(flag,'show_class'))
   
    
    for i=iter_num(1):iter_num(2)
        size_stack=align2d(i,1).stack_size;
        file_path=align2d(i,1).file_path;
 
        
        if (strcmp(num_of_parts,'all'))
            numb=align2d(i,1).avg.num_array(class_nr);
         else
            numb=num_of_parts;
         end;
        if (numb > align2d(i,1).avg.num_array(class_nr))
            num=align2d(i,1).avg.num_array(class_nr);
        end;
         
        zz=1;
        st_tmp=zeros(size_stack(1),size_stack(2),numb);
        cccs=zeros(numb,1);
        for ii=1:size(align2d(i,:),2)
            ii
            if (align2d(i,ii).angleclass.proj_nr==class_nr)
                tmp=tom_emreadc([align2d(i,1).filename],'subregion',[1 1 ii],[size_stack(1)-1 size_stack(2)-1 0]);
                st_tmp(:,:,zz)=tmp.Value;
                cccs(zz)=align2d(1,ii).angleclass.ccc_max;
                zz=zz+1;
            end;
            if (zz==numb)
                return;
            end;
        
        end;
   end;
end;










