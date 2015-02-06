function [vol stack]=tom_av2_xmipp_doc2rec_use_mat(doc_filename)
%tom_av2_xmipp_doc2rec_use_mat calculates 3d rec
%
%   [vol stack]=tom_av2_xmipp_doc2rec_use_mat(doc_filename)
%
%   calculates 3d volume and angular classes   
%
%PARAMETERS
%
%  INPUT
%   doc_filename  name of the xmipp *.doc file from ml3d,ProJmatch
%
%  OUTPUT
%  vol      3d reconstruction  
%  stack    angular classes
%
%
%EXAMPLE
%   
% [vol stack]=tom_av2_xmipp_doc2rec_use_mat('Iter_24_current_angles.doc');
%
%REFERENCES
%
%   xmipp
%
%SEE ALSO
%   tom_xmippsellread, tom_xmippdocread
%
%   created by FB 
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
part_stack_alg=zeros(size(im_tmp.Value,1),size(im_tmp.Value,2),size(st,2));
proj_ang2=zeros(3,max([st(i).ref]));

fprintf('%s \n', ['...done! ' ]); 


if (exist('f_alignment','var') && isempty(f_alignment)==0 )
    load(f_alignment);
end;

fprintf('%s ', ['Building Classes: ' ]); 
for i=1:size(st,1)
  
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
            class(i)=st(i).ref;
            part_stack_alg(:,:,i)=im_tmp_alg;
        else
            im.Value=tom_mirror(im.Value,'x');
            im.Value=im.Value;
            im_tmp_alg=tom_rotate(tom_shift(im.Value,[-st(i).xoff st(i).yoff]),tmp_angles(i,3));
            stack2(:,:,st(i).ref)=stack2(:,:,st(i).ref)+im_tmp_alg;
            count2(st(i).ref)=count2(st(i).ref)+1;
            [aa tttemp]=tom_eulerconvert_xmipp(st(i).rot, st(i).tilt, 0);
            proj_ang2(:,st(i).ref)=tttemp;
            class(i)=st(i).ref;
            part_stack_alg(:,:,i)=im_tmp_alg;
        end;
        if (mod(i,200)==0)
            fprintf('%s','.');
        end;
    
end;
fprintf('%s \n', ['...done! ' ]); 


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