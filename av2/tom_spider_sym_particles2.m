function [new_stack st]=tom_spider_sym_particles2(partstack,part_data,proj_angles,sym_opt,sym_opt_flag,clas)

flag='notest';

angles(:,1)=proj_angles(:,4);
angles(:,2)=proj_angles(:,5);

if (strcmp(flag,'test'))
   part_refs=1:799; 
else
    part_refs=part_data(:,6);
end

flip=part_data(:,17);

new_stack=zeros(size(partstack,1),size(partstack,2),size(partstack,3).*2);
%new_stack(:,:,1:size(partstack,3))=partstack;


%transform sym-opt in spider-rotion sceme zxz 2 zyz yeees transport!!
if (strcmp(sym_opt_flag,'zxz'))
    
    %[sym_mat sym_opt]=tom_eulerconvert_xmipp(sym_opt(1),sym_opt(3),sym_opt(2));
end;

ref_nr=zeros(size(partstack,3),1);
all_dists=zeros(size(partstack,3),1);
euler_ang_zxz=zeros(size(partstack,3),3);
euler_ang_zyz=zeros(size(partstack,3),3);

for i=1:size(partstack,3)
    ref_nr(i)=part_refs(i);
    all_dists(i)=0;
end;


for i=1:size(angles,1)
    [xx,euler_ang_zxz(i,:)] = tom_eulerconvert_xmipp(angles(i,2),angles(i,1),0);
    euler_ang_zyz(i,:)=[angles(i,:) 0];
    [xx,tangles] = tom_eulerconvert_xmipp(angles(i,2),angles(i,1),0);
    tproj_angles(i,:) = tangles';
end;


zz=size(partstack,3)+1;

%figure;

for i=1:size(partstack,3)
    if (flip(i)==1)
        new_stack(:,:,i)=partstack(:,:,i);
    else
        new_stack(:,:,i)=tom_mirror(partstack(:,:,i),'x');
        new_stack(:,:,i)=tom_shift(new_stack(:,:,i),[-1 1]);
    end;
end;


for i=1:size(partstack,3)
    %transform in 3d symopt in one in-plane and psi theta check backproj
    %euler
    inplane_opt=180; %to be calculated  ...check tom_backproj3d_euler
    proj_opt=[0 180]; %to be calculated  ...check tom_backproj3d_euler 
    
    
    %apply in-plane operation
    if (flip(i)==1)
        part_rot=tom_rotate(partstack(:,:,i),inplane_opt);
    else
        part_rot=tom_rotate(tom_mirror(partstack(:,:,i),'x'),inplane_opt);
        part_rot=tom_shift(part_rot,[-1 1]);
    end;
    %filp to compensate for southern sphere phi=180
    part_rot_filp=tom_mirror(part_rot,'y');
    part_rot_filp=tom_shift(part_rot_filp,[-1 1]);
    
    %use psi and theta 
    if (angles(part_refs(i),2)-proj_opt(2) < 0)
        newAng=angles(part_refs(i),:)-proj_opt;
    end;
    
    if (angles(part_refs(i),2)-proj_opt(2) > 0)
        newAng=angles(part_refs(i),:)-proj_opt;
        newAng=[0 360]-newAng;
    end;
    
     % newAng=abs(newAng);
    [pointidx, pointcoords, distance] = tom_nearestpoint(abs(newAng),angles);
    new_stack(:,:,zz)=part_rot_filp;
    ref_nr(zz)=pointidx;
    all_dists(zz)=distance;
    [xx,euler_ang_zxz(zz,:)] = tom_eulerconvert_xmipp(pointcoords(1),pointcoords(2),0);
    euler_ang_zyz(zz,:)=[pointcoords 0];
    
    zz=zz+1;
    
%     %for db
%         %apply in-plane operation
%         cl_rot=tom_rotate(clas(:,:,part_refs(i)),inplane_opt);
%         %filp to compensate for southern sphere phi=180
%         cl_rot_filp=tom_mirror(cl_rot,'y');
%         cl_rot_filp=tom_shift(cl_rot_filp,[-1 1]);
%     %end
%     if (flip(i)==1)
%          subplot(2,3,1); tom_imagesc(tom_filter(partstack(:,:,i),3)); title('org part');
%      else
%         subplot(2,3,1); tom_imagesc(tom_filter(tom_mirror(partstack(:,:,i),'y'),3)); title('org part');
%      end;
%     subplot(2,3,2); tom_imagesc(clas(:,:,part_refs(i))); title('org class');
%     subplot(2,3,4); tom_imagesc(tom_filter(part_rot_filp,3)); title('rot flip part');
%     subplot(2,3,5); tom_imagesc(clas(:,:,pointidx)); title('new trans class');
%     subplot(2,3,6); tom_imagesc(cl_rot_filp); title('org class rot flip');
%     % subplot(6,2,3); tom_imagesc(clas(:,:,pointidx)); title('new trans class');
%     %subplot(3,1,3); tom_imagesc(tom_norm(partstack(:,:,pointidx),'mean0+1std')-tom_norm(part_rot_filp,'mean0+1std'));
%     set(gcf,'Name',num2str(newAng));
%     drawnow;  
%     disp(['angle is: ' num2str(angles(i,1)) ' ' num2str(angles(i,2)) ' transformed ' num2str(abs(newAng(1))) ' '   num2str(abs(newAng(2)))  ' corresponding: ' num2str(pointcoords(1)) ' ' num2str(pointcoords(2)) ' flip ' num2str(flip(i))]  );
%     %cc=tom_corr(partstack(:,:,pointidx),tom_shift(part_rot_filp,[0 0]),'norm'); [a b]=tom_peak(cc)
%     disp(' ');
%     
    disp(num2str(zz));
end;


st.ref_nr=ref_nr;
st.all_dists=all_dists;
st.euler_ang_zxz=euler_ang_zxz;
st.euler_ang_zyz=euler_ang_zyz;
st.tproj_angles=tproj_angles;

