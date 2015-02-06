function [vol_var,var_stack,vol,stack,big_stack,indx]=tom_xmippvariancemap3d(doc_filename,f_ind,f_alignment,sym)
%tom_xmippvariancemap3d calculates a 3d variance Map 
%
%   [vol_var,var_stack,vol,stack]=tom_xmippvariancemap3d(doc_filename)
%
%   calculates 2d variance volume by calculation class variances
%   and weighted backprojection
%
%PARAMETERS
%
%  INPUT
%   doc_filename  name of the xmipp *.doc file from ml3d,ProJ match
%   ind_out       filename of parts per class index
%   f_alignament  alignment filename  
%   symm          (C1) Cn
%
%  OUTPUT
%   vol_var           variance volume
%   var_stack         class variances
%   vol               mean volume            
%   stack             mean classes
%   big_stack         proj,mean and variance Stack    
%
%EXAMPLE
%    [vol_var,var_stack,vol,stack,big_stack]=tom_xmippvariancemap3d('26S_out__it00012doc');
%
%REFERENCES
%
%   xmipp
%
%SEE ALSO
%   tom_xmippsellread, tom_xmippdocread
%
%   created by FB (feat) Heinz Schenk 27/07/07
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



if (nargin<4)
    sym='C1';
    numc=1;
    sym_angs(1)=0;
end;


if (strcmp(sym,'C1')==0)
    numc=str2double(sym(2));
    sym_angs(1)=0;
    for i=2:numc;
        sym_angs(i)=360./numc;
    end;
end;

fprintf('\n%s ', ['Reading  ' doc_filename ':']);  
st=tom_xmippdocread(doc_filename);
fprintf('%s \n', ['...done!']); 
disp('');


fprintf('%s ', ['Converting Angles: ' ]); 
for i=1:size(st,1)
    tmp(i)=st(i).ref;
    [xx,angles] = tom_eulerconvert_xmipp(st(i).rot, st(i).tilt, st(i).psi);
    proj_angles(i,:) = angles';
    tmp_angles(i,:)=[st(i).rot st(i).tilt st(i).psi];
    if (mod(i,2000)==0)    
        fprintf('%s','.');
    end;
end;
fprintf('%s \n','done!');



fprintf('%s ', ['Allocating Memory: ' ]); 
num_of_ref=max(tmp);
try
    im_tmp=tom_spiderread(['./' st(1).name]);
     path_flag='rel_man';
catch
    im_tmp=tom_spiderread([st(1).name]);
    path_flag='abs_man';
end


stack=zeros(size(im_tmp.Value,1),size(im_tmp.Value,2),num_of_ref);
stack2=zeros(size(im_tmp.Value,1),size(im_tmp.Value,2),num_of_ref);
count=zeros(num_of_ref);
count2=zeros(num_of_ref);
class=zeros(size(st,2),1);
feat_stack=zeros(size(st,2),size(im_tmp.Value,1).*size(im_tmp.Value,1));
part_stack_alg=zeros(size(im_tmp.Value,1),size(im_tmp.Value,2),size(st,2));

%indx_p=zeros(num_of_ref,1);
for i=1:num_of_ref
    indx(i).filename{1}='';
%     indx(i).rot(1)='';
%     indx(i).shift(1,1)='';
%     indx(i).shift(1,2)='';
end;

fprintf('%s \n', ['...done! ' ]); 


if (exist('f_alignment','var') && isempty(f_alignment)==0 )
    load(f_alignment);
end;

max_sh=0.10.*size(im_tmp.Value,1);
sort_out_count=0;
fprintf('%s ', ['Building Classes: ' ]); 
for i=1:size(st,1)
    if (abs(st(i).xoff) < max_sh && abs(st(i).yoff) < max_sh )
        
        if (strcmp(path_flag,'rel_man'))
            tmp_name=['./' st(i).name];
        else
            tmp_name=[st(i).name];
        end;
        try
            im=tom_spiderread(tmp_name);
        catch
            disp(['particle not found ' tmp_name]);
            continue;
        end;
        
        im.Value=tom_xraycorrect(im.Value);
        im.Value=tom_norm(im.Value,'mean0+1std');
        
        if (st(i).flip==0)
            im_tmp_alg=tom_rotate(tom_shift(im.Value,[st(i).xoff st(i).yoff]),tmp_angles(i,3));
            stack(:,:,st(i).ref)=stack(:,:,st(i).ref)+im_tmp_alg;
            count(st(i).ref)=count(st(i).ref)+1;
            [aa tttemp]=tom_eulerconvert_xmipp(st(i).rot, st(i).tilt, 0);
            proj_ang2(:,st(i).ref)=tttemp;
            %feat_stack(i,:)=reshape(tom_norm(im_tmp_alg,'mean0+1std',mask).*mask,size(im_tmp.Value,1).*size(im_tmp.Value,1),1);
            feat_stack(i,:)=reshape(im_tmp_alg,size(im_tmp.Value,1).*size(im_tmp.Value,1),1);
            class(i)=st(i).ref;
            
            if (isempty(indx(st(i).ref).filename{1}))
                tmp_indd=length(indx(st(i).ref).filename);
            else
                tmp_indd=length(indx(st(i).ref).filename)+1;
            end;
            indx(st(i).ref).filename{tmp_indd}=tmp_name;
            indx(st(i).ref).rot(tmp_indd)=tmp_angles(i,3);
            indx(st(i).ref).shift(tmp_indd,:)=[st(i).xoff st(i).yoff];
            indx(st(i).ref).flip(tmp_indd)=st(i).flip;
            indx(st(i).ref).sel(tmp_indd)=1;
            if (exist('f_alignment','var') && isempty(f_alignment)==0)
                indx(st(i).ref).align{tmp_indd}=align2d(1,i);
            end;
            part_stack_alg(:,:,i)=im_tmp_alg;
        else
            im.Value=tom_mirror(im.Value,'x');
            im.Value=im.Value;
            im_tmp_alg=tom_rotate(tom_shift(im.Value,[-st(i).xoff st(i).yoff]),tmp_angles(i,3));
            stack2(:,:,st(i).ref)=stack2(:,:,st(i).ref)+im_tmp_alg;
            count2(st(i).ref)=count2(st(i).ref)+1;
            [aa tttemp]=tom_eulerconvert_xmipp(st(i).rot, st(i).tilt, 0);
            proj_ang2(:,st(i).ref)=tttemp;
            %feat_stack(i,:)=reshape(tom_norm(im_tmp_alg,'mean0+1std',mask).*mask,size(im_tmp.Value,1).*size(im_tmp.Value,1),1);
            feat_stack(i,:)=reshape(im_tmp_alg,size(im_tmp.Value,1).*size(im_tmp.Value,1),1);
            class(i)=st(i).ref;
            
            if (isempty(indx(st(i).ref).filename{1}))
                tmp_indd=length(indx(st(i).ref).filename);
            else
                tmp_indd=length(indx(st(i).ref).filename)+1;
            end;
            indx(st(i).ref).filename{tmp_indd}=tmp_name;
            indx(st(i).ref).rot(tmp_indd)=tmp_angles(i,3);
            indx(st(i).ref).shift(tmp_indd,:)=[-st(i).xoff st(i).yoff];
            indx(st(i).ref).flip(tmp_indd)=st(i).flip;
            indx(st(i).ref).sel(tmp_indd)=1;
            if (exist('f_alignment','var') && isempty(f_alignment)==0)
                indx(st(i).ref).align{tmp_indd}=align2d(1,i);
            end;
            part_stack_alg(:,:,i)=im_tmp_alg;
        end;
        if (mod(i,200)==0)
            fprintf('%s','.');
        end;
    else
        sort_out_count=sort_out_count+1;
    end;
end;
fprintf('%s \n', ['...done! ' ]); 

disp([num2str(sort_out_count) ' particles sorted out']);

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

mask=tom_spheremask(ones(size(stack,1),size(stack,2)),round(size(stack,1)./2));



fprintf('%s ', ['Building Varianz Classes: ' ]); 
var_stack=zeros(size(im_tmp.Value,1),size(im_tmp.Value,1),max(class));
%build up Variance stack
parfor i=1:max(class)
%for i=1:max(class)
 %   var_stack(:,:,i)=tom_norm(reshape(tom_calc_variance(feat_stack,class,i),size(im_tmp.Value,1),size(im_tmp.Value,1))+100,'mean0+1std');
    tmp_cl_stack=part_stack_alg(:,:,find(class==i));
    if (isempty(tmp_cl_stack) || size(tmp_cl_stack,3)==1 )
        var_stack(:,:,i)=zeros(size(mask,1),size(mask,2));
    else
        var_stack(:,:,i)=tom_calc_variance2d(tmp_cl_stack,mask,'mean0+1std');
    end;
    if (mod(i,50)==0)    
        fprintf('%s','.');
        drawnow;
    end;
end;
fprintf('%s \n', ['...done! ' ]); 


stack=stack+stack2;
%adapt angles
zz=1;
tmpp=zeros(3,size(stack,3));
for i=1:size(stack,3)
    if (std2(stack(:,:,i) )~=0 )
        for ii=1:numc
            tmpp(:,zz)=[(proj_ang2(1,i)+sym_angs(ii)) proj_ang2(2,i) proj_ang2(3,i)];
            zz=zz+1;
        end;
    end;
end;
proj_angles_clean=tmpp(:,1:(zz-1));


fprintf('%s ', ['Reconstructing: ' ]); 
%backproject image stack

vol=zeros(size(stack,1),size(stack,1),size(stack,1),'single');
THICK=size(vol,3);

for i=1:size(proj_ang2,2)
   if (std2(stack(:,:,i) )~=0 )
        for ii=1:numc
            w_func=tom_calc_weight_functionc([size(vol,1) size(vol,2)],proj_angles_clean',THICK,[(proj_ang2(1,i)+sym_angs(ii)) proj_ang2(2,i) proj_ang2(3,i)]);
            w_proj=tom_apply_weight_function(stack(:,:,i),w_func);
           
            tom_backproj3d_euler(vol,w_proj,(proj_ang2(1,i)+sym_angs(ii)),proj_ang2(2,i),proj_ang2(3,i),[0 0 0]);
        end;
        if (mod(i,50)==0)
            fprintf('%s','.');
        end;
    end;
    
end;
fprintf('%s \n', ['...done! ' ]); 


%adapt angles
zz=1;
tmpp=zeros(3,size(stack,3));
for i=1:size(var_stack,3)
    if (std2(var_stack(:,:,i) )~=0 )
        for ii=1:numc
            tmpp(:,zz)=[(proj_ang2(1,i)+sym_angs(ii)) proj_ang2(2,i) proj_ang2(3,i)];
            zz=zz+1;
        end;
    end;
end;
proj_angles_clean=tmpp(:,1:(zz-1));



fprintf('%s ', ['Reconstructing Variance Map: ' ]); 
%backproject variance stack
vol_var=zeros(size(stack,1),size(stack,1),size(stack,1),'single');
for i=1:size(proj_ang2,2)
    if (std2(var_stack(:,:,i) )~=0 )
        for ii=1:numc
            w_func=tom_calc_weight_functionc([size(vol_var,1) size(vol_var,2)],proj_angles_clean',THICK,[(proj_ang2(1,i)+sym_angs(ii)) proj_ang2(2,i) proj_ang2(3,i)]);
            w_proj=tom_apply_weight_function(var_stack(:,:,i),w_func);
           % w_proj=var_stack(:,:,i);
            tom_backproj3d_euler(vol_var,w_proj,(proj_ang2(1,i)+sym_angs(ii)),proj_ang2(2,i),proj_ang2(3,i),[0 0 0]);
        end;
        
        if (mod(i,50)==0)
            fprintf('%s','.');
        end;
    end;
end;
fprintf('%s \n', ['...done! ' ]); 


fprintf('%s ', ['Reprojecting : ' ]); 
proj_stack=zeros(size(im_tmp.Value,1),size(im_tmp.Value,1),max(class));
for i=1:size(proj_ang2,2)
    vol_rot=tom_rotate(vol,[proj_ang2(1,i) proj_ang2(2,i) proj_ang2(3,i)]);
    proj_stack(:,:,i)=sum(vol_rot,3);
    if (mod(i,50)==0)    
        fprintf('%s','.');
    end;
end;

fprintf('%s \n', ['...done! ' ]); 

fprintf('%s ', ['Combining Stack : ' ]); 



big_stack=zeros(size(im_tmp.Value,1),size(im_tmp.Value,1),max(class).*3);
zz=1;
for i=1:3:((size(proj_ang2,2)*3)-3)
    big_stack(:,:,i)=tom_norm(proj_stack(:,:,zz),'mean0+1std');
    big_stack(:,:,i+1)=tom_norm(tom_filter(stack(:,:,zz),1),'mean0+1std');
    big_stack(:,:,i+2)=tom_norm(tom_filter(var_stack(:,:,zz),1),'mean0+1std',mask).*mask;
    zz=zz+1;
    if (mod(i,50)==0)
        fprintf('%s','.');
    end;
end;
fprintf('%s \n', ['...done! ' ]); 

if (exist('f_ind','var')  && isempty(f_ind)==0)
    try
        save(f_ind,'indx','big_stack');
        clear f_ind;
    catch
        disp('Cannot write comp');
    end;
end;

clear f_ind;




