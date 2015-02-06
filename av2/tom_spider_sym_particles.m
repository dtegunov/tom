function tom_spider_sym_particles(partstack,proj_angles,sym_opt,sym_opt_flag)

angles(:,1)=proj_angles(:,4);
angles(:,2)=proj_angles(:,5);


%transform sym-opt in spider-rotion sceme zxz 2 zyz yeees transport!!
if (strcmp(sym_opt_flag,'zxz'))
    %sym_opt=tom_eulerconvert_xmipp(sym_opt);
end;

figure; 

for i=622:size(partstack,3)
    %transform in 3d symopt in one in-plane and psi theta chech backproj
    %euler
    inplane_opt=180; %to be calculated  ...check tom_backproj3d_euler
    proj_opt=[0 180]; %to be calculated  ...check tom_backproj3d_euler 
    
    
    %apply in-plane operation
    part_rot=tom_rotate(partstack(:,:,i),inplane_opt);
    %filp to compensate for southern sphere phi=180
    part_rot_filp=tom_mirror(part_rot,'y');
    part_rot_filp=tom_shift(part_rot_filp,[-1 1]);
    
    %use psi and theta 
    if (angles(i,2)-proj_opt(2) < 0)
        newAng=angles(i,:)-proj_opt;
    end;
    
    if (angles(i,2)-proj_opt(2) > 0)
        newAng=angles(i,:)-proj_opt;
        newAng=[0 360]-newAng;
    end;
    
%      if (angles(i,2)==0 )
%         newAng=angles(i,:)-proj_opt;
%         %newAng=[0 360]-newAng;
%     end;
    
   % newAng=abs(newAng);
    [pointidx, pointcoords, distance] = tom_nearestpoint(abs(newAng),angles);
    
    subplot(3,1,1); tom_imagesc(part_rot_filp);
    subplot(3,1,2); tom_imagesc(partstack(:,:,pointidx));
    subplot(3,1,3); tom_imagesc(tom_norm(partstack(:,:,pointidx),'mean0+1std')-tom_norm(part_rot_filp,'mean0+1std'));
    set(gcf,'Name',num2str(newAng));
    drawnow;  
    disp(['angle is: ' num2str(angles(i,1)) ' ' num2str(angles(i,2)) ' transformed ' num2str(abs(newAng(1))) ' '   num2str(abs(newAng(2)))  ' corresponding: ' num2str(pointcoords(1)) ' ' num2str(pointcoords(2))]);
    %cc=tom_corr(partstack(:,:,pointidx),tom_shift(part_rot_filp,[0 0]),'norm'); [a b]=tom_peak(cc)
    disp(' ');
    
end;

