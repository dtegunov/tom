function [classes vol classes_var vol_var]=tom_spidervariancemap3d(alignedStackName,paramfileName,proj_ang)

disp('Reading data ...');

stack_alg=tom_readspiderstack(alignedStackName);
stack_alg=stack_alg.Value;

data=importdata(paramfileName);
data=data.data;

proj_data=importdata(proj_ang);
proj_data=proj_data.data;

disp('Allocating memory ...')
ref_nr=data(:,6);
psi=data(:,4);
theta=data(:,5);
rot=data(:,8);
flip=data(:,17);
classes=zeros(size(stack_alg,1),size(stack_alg,2),max(ref_nr));
part_per_class=zeros(max(ref_nr),1);
proj_angles=zeros(size(classes,3),3);
psi=proj_data(:,4);
theta=proj_data(:,5);


fprintf('%s ', 'Building classes: ' ); 
flip_stack_alg=zeros(size(stack_alg));


for i=1:size(stack_alg,3)
    if (flip(i)==1)
        classes(:,:,ref_nr(i))=classes(:,:,ref_nr(i))+stack_alg(:,:,i);
        flip_stack_alg(:,:,i)=stack_alg(:,:,i);
     else
       classes(:,:,ref_nr(i))=classes(:,:,ref_nr(i))+tom_mirror(stack_alg(:,:,i),'x');
       flip_stack_alg(:,:,i)=tom_mirror(stack_alg(:,:,i),'x');
     end;
     part_per_class(ref_nr(i))=part_per_class(ref_nr(i))+1;
    if (mod(i,1000)==0)    
        fprintf('%s','.');
    end;
   
end;
fprintf('%s \n', '...done! '); 


fprintf('%s ', 'Norming classes: ' ); 
for i=1:size(classes,3)
   if (part_per_class(i)==0)
       part_per_class(i)=1;
    end;
    classes(:,:,i)=tom_norm(classes(:,:,i)./part_per_class(i),'mean0+1std');
    if (mod(i,20)==0)    
        fprintf('%s','.');
    end;
   [xx,angles] = tom_eulerconvert_xmipp(theta(i),psi(i),0);
   proj_angles(i,:) = angles';
end;
fprintf('%s \n', '...done! '); 



fprintf('%s ', ['Reconstructing volume: ' ]); 
%backproject image stack
vol=zeros(size(stack_alg,1),size(stack_alg,1),size(stack_alg,1),'single');
for i=1:size(proj_angles,1)
    w_func=tom_calc_weight_functionc([size(vol,1) size(vol,2)],proj_angles,80,[proj_angles(i,1) proj_angles(i,2) proj_angles(i,3)]);
    w_proj=tom_apply_weight_function(classes(:,:,i),w_func);
    tom_backproj3d_euler(vol,w_proj,proj_angles(i,1),proj_angles(i,2),proj_angles(i,3),[0 0 0]);
    if (mod(i,50)==0)    
        fprintf('%s','.');
    end;
end;
figure; tom_dspcub(vol,1);
figure; tom_dspcub(vol,2);
fprintf('%s \n', ['...done! ' ]); 


fprintf('%s ', 'building var classes: ' );  
classes_var=zeros(size(stack_alg,1),size(stack_alg,2),max(ref_nr));
for i=1:size(classes,3)
    part_tmp=flip_stack_alg(:,:,find(ref_nr==i));
    classes_var(:,:,i)=tom_calc_variance_stack(part_tmp);
    if (mod(i,50)==0)    
         fprintf('%s','.');
    end;
end;
fprintf('%s \n', '...done! '); 



fprintf('%s ', 'Norming classes: ' ); 
for i=1:size(classes,3)
   if (part_per_class(i)==0)
       part_per_class(i)=1;
    end;
    classes_var(:,:,i)=classes_var(:,:,i)./part_per_class(i);
    if (mod(i,20)==0)    
        fprintf('%s','.');
    end;
end;
fprintf('%s \n', '...done! '); 


fprintf('%s ', ['Reconstructing var volume: ' ]); 
%backproject image stack
vol_var=zeros(size(stack_alg,1),size(stack_alg,1),size(stack_alg,1),'single');
for i=1:size(proj_angles,1)
    w_func=tom_calc_weight_functionc([size(vol,1) size(vol,2)],proj_angles,80,[proj_angles(i,1) proj_angles(i,2) proj_angles(i,3)]);
    w_proj=tom_apply_weight_function(classes_var(:,:,i),w_func);
    tom_backproj3d_euler(vol_var,w_proj,proj_angles(i,1),proj_angles(i,2),proj_angles(i,3),[0 0 0]);
    if (mod(i,50)==0)    
        fprintf('%s','.');
    end;
end;
fprintf('%s \n', ['...done! ' ]); 












