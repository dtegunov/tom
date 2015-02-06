function [stack,st,vol]=tom_xmippdoc2statstack(filename)

st=tom_xmippdocread(filename);


for i=1:size(st,2)
    tmp(i)=st(i).ref;
    [xx,angles] = tom_eulerconvert_xmipp(st(i).rot, st(i).tilt, st(i).psi);
    proj_angles(i,:) = angles';
    tmp_angles(i,:)=[st(i).rot st(i).tilt st(i).psi];
end;

num_of_ref=max(tmp);
im_tmp=tom_spiderread(['.' st(1).name]);


stack=zeros(size(im_tmp.Value,1),size(im_tmp.Value,2),num_of_ref);
stack2=zeros(size(im_tmp.Value,1),size(im_tmp.Value,2),num_of_ref);

count=zeros(num_of_ref);
count2=zeros(num_of_ref);

for i=1:size(st,2)
    im=tom_spiderread(['.' st(i).name]);
    %stack(:,:,st(i).ref)=stack(:,:,st(i).ref)+tom_shift(tom_rotate(im.Value,st(i).psi),[st(i).xoff st(i).yoff]);
    %st(i).flip=0;
    
    if (st(i).ref==1)
        disp('tt');
    end;
    
    if (st(i).flip==0)
        im_tmp_alg=tom_rotate(tom_shift(im.Value,[st(i).xoff st(i).yoff]),tmp_angles(i,3));
        stack(:,:,st(i).ref)=stack(:,:,st(i).ref)+im_tmp_alg;
        count(st(i).ref)=count(st(i).ref)+1;
        [aa tttemp]=tom_eulerconvert_xmipp(st(i).rot, st(i).tilt, st(i).psi);
        proj_ang2(:,st(i).ref)=tttemp;
    else
        im.Value=tom_mirror(im.Value,'x');
        im.Value=im.Value;
        im_tmp_alg=tom_rotate(tom_shift(im.Value,[-st(i).xoff st(i).yoff]),tmp_angles(i,3));
        stack2(:,:,st(i).ref)=stack2(:,:,st(i).ref)+im_tmp_alg;
        count2(st(i).ref)=count2(st(i).ref)+1;
        [aa tttemp]=tom_eulerconvert_xmipp(st(i).rot, st(i).tilt, st(i).psi);
        proj_ang2(:,st(i).ref)=tttemp;
    end;
   
    disp(num2str(i));
end;

%norm it norman
for i=1:size(stack,3)
    if (count(i)~=0)
        stack(:,:,i)=tom_norm(stack(:,:,i)./count(i)+100,'phase');
    end;
    
end;

for i=1:size(stack2,3)
    if (count2(i)~=0) 
        stack2(:,:,i)=tom_norm(stack2(:,:,i)./count2(i)+100,'phase');
    end;
end;

%stack=stack+stack2;

for i=1:size(stack,3)
    stack(:,:,i)=tom_norm(stack(:,:,i)+100,'phase');
    
end;
%stack = '';
%load wss


vol=zeros(80,80,80,'single');

%vol=tom_bin(vol,1);




for i=1:size(stack,3)
    %if (mod(i,2)==1)
    %proj=stack(:,:,i);
    % proj=tom_bin(proj,1);
    %w_func=tom_calc_weight_functionc(round(size(im_tmp.Value)./2),tmp_angles,160,[0 proj_angles(1,i) proj_angles(2,i)]);
    % w_proj=tom_apply_weight_function(proj,w_func);
    %w_proj=proj;
    %if (i==438)
    %    disp('emd');
    %end;

    
    %w_func=tom_calc_weight_functionc([size(vol,1) size(vol,2)],proj_ang2',80,[proj_ang2(1,i) proj_ang2(2,i) proj_ang2(3,i)]);
    %w_proj=tom_apply_weight_function(stack(:,:,i),w_func);
     
    w_proj=single(stack(:,:,i));
    
    
    %tom_backproj3d_euler(vol,w_proj,proj_ang2(1,i),proj_ang2(3,i),proj_ang2(2,i),[0 0 0]);
    %tom_backproj3d_euler(vol,w_proj,0,proj_ang2(3,i),proj_ang2(2,i),[0 0 0]);
    tom_backproj3d(vol,w_proj,proj_ang2(2,i),proj_ang2(3,i),[0 0 0]);
    disp(i);
    
end;







return;



for i=1:size(st,2)
    %if (mod(i,2)==1)
    %proj=stack(:,:,i);
    % proj=tom_bin(proj,1);
    %w_func=tom_calc_weight_functionc(round(size(im_tmp.Value)./2),tmp_angles,160,[0 proj_angles(1,i) proj_angles(2,i)]);
    % w_proj=tom_apply_weight_function(proj,w_func);
    %w_proj=proj;
    %if (i==438)
    %    disp('emd');
    %end;

    im=tom_spiderread(['.' st(i).name]);
    if st(i).flip==1
%        continue;
        im.Value=tom_mirror(im.Value,'x');
        st(i).xoff = -st(i).xoff;
    end
    w_proj = tom_shift(im.Value,[st(i).xoff st(i).yoff]);
    w_func=tom_calc_weight_functionc(round(size(im_tmp.Value)),proj_angles,48,[proj_angles(i,1) proj_angles(i,2) proj_angles(i,3)]);
    w_proj=tom_apply_weight_function(w_proj,w_func);
    
    tom_backproj3d_euler(vol,w_proj,proj_angles(i,1),proj_angles(i,2),proj_angles(i,3),[0 0 0]);
    disp(i);
    
end;



disp('end');