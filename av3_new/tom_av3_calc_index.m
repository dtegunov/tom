function tom_av3_calc_index(stack_in,path_out,ilevel, iknoten, part,name,flag)


%tom_av3_calc_index(im.Value,'data',1,2,[1:72],[],'classify');


if strcmp(flag,'2_halfs')==1 
    N = ceil(length(part) / 2);
    part_left = part(1:N);
    part_right = part(1+N:end);
end;

if strcmp(flag,'classify')==1 
    idx=tom_pca_and_k_means(stack_in(:,:,part),2,5,1,0);
    part_left=part(idx==1);
    part_right=part(idx==2);
end;





if (isempty(part_left))
else
    name(length(name)+1)=0;
    disp(['members: ' num2str(part_left) ' name: ' num2str(name) ]);
    mean_tmp=zeros(size(stack_in,1),size(stack_in,2));
    mean_tmp=tom_emheader(mean_tmp);
    for i=1:length(part_left)
        mean_tmp.Value=mean_tmp.Value+stack_in(:,:,part_left(i));
    end;
    
    mean_tmp.Header.Comment=char(strrep(num2str(part_left),' ',''));
    t_name=strrep(num2str(name),' ','');
    tom_emwrite([path_out '/idx_' t_name '.em'],mean_tmp);    
    
    if (length(part_left) > 1)
        tom_av3_calc_index(stack_in,path_out,ilevel+1, iknoten*2+0, part_left,name,flag)
    end;
end;



if (isempty(part_right))
else
    name(length(name))=1;
    disp(['members ' num2str(part_right) ' name: ' num2str(name) ]);
    if (strcmp(strrep(num2str(name),' ','') ,'01'))
        disp('');
    end;
    mean_tmp=zeros(size(stack_in,1),size(stack_in,2));  
    mean_tmp=tom_emheader(mean_tmp);
    for i=1:length(part_right)
        mean_tmp.Value=mean_tmp.Value+stack_in(:,:,part_right(i));
    end;
    mean_tmp.Header.Comment=char(strrep(num2str(part_right),' ',''));
    t_name=strrep(num2str(name),' ','');
    tom_emwrite([path_out '/idx_' t_name '.em'],mean_tmp);   
    
    if (length(part_right) > 1)
        tom_av3_calc_index(stack_in,path_out,ilevel+1, iknoten*2+1, part_right,name,flag)
    end;
end;






