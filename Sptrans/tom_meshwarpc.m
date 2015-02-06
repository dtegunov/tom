function [image_out,mesh1_out,mesh2_out] =tom_meshwarpc(image_in,mesh1,mesh2)
%TOM_MESHWARPC performs a warping by calling a c-function
%
%   [image_out,mesh1_out,mesh2_out] =tom_meshwarpc(image_in,mesh1,mesh2)
%
%   Performs worping using the algorithm of George Wolberg which is implemented 
%   in the C-Function(mashworp.c). The matlab fuctions transforms the input
%   meshes (mesh1 mesh2) in the format wolber algorithm useses and resolves
%   conflicts caused by the movment of the markers.
%
%PARAMETERS
%
%  INPUT
%   image_in            input image
%   mesh1               mesh with the start markers
%   mesh2               mesh with the moved markers
%  
%  OUTPUT
%   image_out           warped image
%   mesh1_out           mesh with the start markers without conflicts
%   mesh2_out           mesh with the moved markers without conflicts
%
%EXAMPLE
%   [imout m1 m2]=tom_meshwarpc(cat,base_points,input_points); figure; tom_imagesc(imout);
%
%   mesh1    mesh2
%    x   y    x   y
%    x   y    x   y 
%    x   y    x   y 
%
%NOTE
%   The input meshes must have the following form
%
%REFERENCES
%
%SEE ALSO
%   TOM_UNBEND, CPSELECT
%
%   created by FB 04/13/04
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

%cast variables
image_in=single(image_in);
mesh1=single(mesh1);
mesh2=single(mesh2);
image_out=single(zeros(size(image_in)));

%begin of the loop to resolve conflicts of x movement
%the first marker causing a confilct is deleted
conflict=1;
conf_z=1;
while (conflict==1)
    
    % calculate the dimensions of the matrix mesh_xy1 and mesh_xy2
    n=size(mesh1,1); 
    dim=floor(sqrt(n));
    rest=n-(dim.*dim);
    if (rest==0)
        dim=dim;
    else
        dim=dim+1;
    end;
    rest=(dim.*dim)-n;
    
    %allocate memory
    mesh_xy1=single(zeros(dim+2,dim+2,2));
    mesh_xy2=single(zeros(dim+2,dim+2,2));
    
    % add markers
    inkr=2.12452; % strange increment to avoid conflicts with other markers
    if (rest~=0);
        for i=1:(rest)
            mesh1(size(mesh1,1)+1,:)=[size(image_in,1)-inkr size(image_in,1)-inkr];
            mesh2(size(mesh2,1)+1,:)=[size(image_in,1)-inkr size(image_in,1)-inkr];
        end;
    end;
    
    %sort the meshes by y of mesh1
    tmp=mesh1;
    tmp(:,3:4)=mesh2;
    tmp=sortrows(tmp,2);
    mesh1=tmp(:,1:2);
    mesh2=tmp(:,3:4);
    
    
    %calculate check distance to find conflicts
    check_dist_old=mesh1-mesh2;
    
    % add check distance to mesh1 to have it in the right order
    mesh1(:,3:4)=check_dist_old;
    % sort every line by x
    for pst=1:dim:n
        tmp=mesh1(pst:pst+dim-1,:);
        tmp=sortrows(tmp,1);
        mesh1(pst:pst+dim-1,:)=tmp;
        tmp=mesh2(pst:pst+dim-1,:);
        tmp=sortrows(tmp,1);
        mesh2(pst:pst+dim-1,:)=tmp;
    end;
    
    %check for conflicts
    check_dist_old=mesh1(:,3:4);
    mesh1=mesh1(:,1:2);
    check_dist_new=mesh1-mesh2;
    index=find(check_dist_old(:,1)-(check_dist_new(:,1)));
    if (isempty(index))
        conflict=0;    
    else
        % log the conflicts
        conflict_table(conf_z,1:2)=mesh1(index(1),:);
        conflict_table(conf_z,3:4)=mesh2(index(1),:);
        
        %delete the first conflict marker by sorting and reshaping the
        %matrix
        mesh1(index(1),:)=[-1 -1];
        mesh2(index(1),:)=[-1 -1];
        tmp1=sortrows(mesh1,2);
        tmp2=sortrows(mesh2,2);
        mesh1=tmp1(2:size(tmp1,1),:);
        mesh2=tmp2(2:size(tmp2,1),:);  
        conf_z=conf_z+1;    
    end;
end;

% Output to matlab
if (conf_z > 1)
disp('  ');
fprintf('deleted %d markers to avoid conflicts\n',(conf_z-1));
fprintf('\n         source                 dest \n');
disp(conflict_table);
else
   disp('no conflicts accured'); 
end;

% put the sorted meshes in a different format (x(n,n) y(n,n)) 
zz=1;
for i=2:dim+1
    for ii=2:dim+1
        mesh_xy1(ii,i,1)=mesh1(zz,1);
        mesh_xy1(ii,i,2)=mesh1(zz,2);
        mesh_xy2(ii,i,1)=mesh2(zz,1);
        mesh_xy2(ii,i,2)=mesh2(zz,2);
        zz=zz+1;
    end;
    
end;

%create fix points 
vy=[(size(image_in,2)-1) (size(image_in,2)-1) (size(image_in,2)-1)];
vx=[size(image_in,1)-10];
for i=1:dim-1
    vy=[vy size(image_in,2)-1];
    vx=[vx (size(image_in,1)-i*(((size(image_in,1)./(dim))-2)))];
end;
vx=sort(vx);

for i=2:dim+1
   vy_var(i-1)=mean(mesh_xy1(2:size(mesh_xy1,1)-1,i,2));
end;

mesh_xy1(size(mesh_xy1,1),:,1)=vy;
mesh_xy2(size(mesh_xy1,1),:,1)=vy;

mesh_xy1(:,size(mesh_xy1,2),2)=vy;
mesh_xy2(:,size(mesh_xy1,2),2)=vy;

mesh_xy1(2:(size(mesh_xy1,2)-1),1,1)=vx;
mesh_xy2(2:(size(mesh_xy1,2)-1),1,1)=vx;

mesh_xy1(2:size(mesh_xy1,1)-1,size(mesh_xy1,1),1)=vx;
mesh_xy2(2:size(mesh_xy1,1)-1,size(mesh_xy1,1),1)=vx;

mesh_xy1(1,2:size(mesh_xy1,2)-1,2)=vy_var;
mesh_xy2(1,2:size(mesh_xy1,2)-1,2)=vy_var;

mesh_xy1(size(mesh_xy1,1),2:size(mesh_xy1,2)-1,2)=vy_var;
mesh_xy2(size(mesh_xy1,1),2:size(mesh_xy1,2)-1,2)=vy_var;


%cast mash to single 
mesh_xy1=single(mesh_xy1);
mesh_xy2=single(mesh_xy2);

%image_in=char(image_in);

% do the actual worping
meshwarp(image_in,mesh_xy1,mesh_xy2,image_out);

mesh1_out=mesh1;
mesh2_out=mesh2;
%meshwarp(image_in,mesh1,mesh2,image_out);