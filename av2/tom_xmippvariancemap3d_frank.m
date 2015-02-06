function [vol_var,var_stack,vol,stack,all_vols,all_vars,var_3d,vol_sum]=tom_xmippvariancemap3d_frank(doc_filename,num_of_models)
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
%   doc_filename  name of the xmipp *.doc file from ml3d
%   
%  
%  OUTPUT
%   vol_var           variance volume
%   var_stack         class variances
%   vol               mean volume            
%   stack             mean classes
%
%EXAMPLE
%    [vol_var,var_stack,vol,stack]=tom_xmippvariancemap3d('26S_out__it00012doc');
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



fprintf('\n%s ', ['Reading  ' doc_filename ':']);  
st_all=tom_xmippdocread(doc_filename);
fprintf('%s \n', ['...done!']); 
disp('');

packages=tom_calc_packages(num_of_models,length(st_all));



for iiii=1:num_of_models

    disp(['Building ' num2str(iiii) ' Volume']);
    
    %st=st_all(packages(iiii,1):packages(iiii,2));
    
    idx=round(length(st_all).*rand(8500,1));
    idx=idx+(idx==0);
    st=st_all(idx);
    
    fprintf('%s ', ['Converting Angles: ' ]);
    for i=1:size(st,2)
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
    %im_tmp=tom_spiderread(['./' st(1).name]);
    im_tmp=tom_spiderread([ st(1).name]);
    stack=zeros(size(im_tmp.Value,1),size(im_tmp.Value,2),num_of_ref);
    stack2=zeros(size(im_tmp.Value,1),size(im_tmp.Value,2),num_of_ref);
    count=zeros(num_of_ref);
    count2=zeros(num_of_ref);
    class=zeros(size(st,2),1);
    feat_stack=zeros(size(st,2),size(im_tmp.Value,1).*size(im_tmp.Value,1));
    fprintf('%s \n', ['...done! ' ]);




    fprintf('%s ', ['Building Classes: ' ]);
    for i=1:size(st,2)
        %    im=tom_spiderread(['./' st(i).name]);
        im=tom_spiderread([ st(i).name]);
        if (st(i).flip==0)
            im_tmp_alg=tom_rotate(tom_shift(im.Value,[st(i).xoff st(i).yoff]),tmp_angles(i,3));
            stack(:,:,st(i).ref)=stack(:,:,st(i).ref)+im_tmp_alg;
            count(st(i).ref)=count(st(i).ref)+1;
            [aa tttemp]=tom_eulerconvert_xmipp(st(i).rot, st(i).tilt, 0);
            proj_ang2(:,st(i).ref)=tttemp;
            feat_stack(i,:)=reshape(im_tmp_alg,size(im_tmp.Value,1).*size(im_tmp.Value,1),1);
            class(i)=st(i).ref;
        else
            im.Value=tom_mirror(im.Value,'x');
            im.Value=im.Value;
            im_tmp_alg=tom_rotate(tom_shift(im.Value,[-st(i).xoff st(i).yoff]),tmp_angles(i,3));
            stack2(:,:,st(i).ref)=stack2(:,:,st(i).ref)+im_tmp_alg;
            count2(st(i).ref)=count2(st(i).ref)+1;
            [aa tttemp]=tom_eulerconvert_xmipp(st(i).rot, st(i).tilt, 0);
            proj_ang2(:,st(i).ref)=tttemp;
            feat_stack(i,:)=reshape(im_tmp_alg,size(im_tmp.Value,1).*size(im_tmp.Value,1),1);
            class(i)=st(i).ref;
        end;
        if (mod(i,2000)==0)
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




    fprintf('%s ', ['Building Varianz Classes: ' ]);
    var_stack=zeros(size(im_tmp.Value,1),size(im_tmp.Value,1));
    %build up Variance stack
    for i=1:max(class)
        var_stack(:,:,i)=tom_norm(reshape(tom_calc_variance(feat_stack,class,i),size(im_tmp.Value,1),size(im_tmp.Value,1))+100,'phase');
        if (mod(i,50)==0)
            fprintf('%s','.');
        end;
    end;
    fprintf('%s \n', ['...done! ' ]);




    fprintf('%s ', ['Reconstructing: ' ]);
    %backproject image stack
    stack=stack+stack2;
    vol=zeros(size(stack,1),size(stack,1),size(stack,1),'single');
    for i=1:size(proj_ang2,2)
        w_func=tom_calc_weight_functionc([size(vol,1) size(vol,2)],proj_ang2',80,[proj_ang2(1,i) proj_ang2(2,i) proj_ang2(3,i)]);
        w_proj=tom_apply_weight_function(stack(:,:,i),w_func);
        tom_backproj3d_euler(vol,w_proj,proj_ang2(1,i),proj_ang2(2,i),proj_ang2(3,i),[0 0 0]);
        if (mod(i,50)==0)
            fprintf('%s','.');
        end;
    end;
    fprintf('%s \n', ['...done! ' ]);


    fprintf('%s ', ['Reconstructing Variance Map: ' ]);
    %backproject variance stack
    vol_var=zeros(size(stack,1),size(stack,1),size(stack,1),'single');
    for i=1:size(proj_ang2,2)
        w_func=tom_calc_weight_functionc([size(vol_var,1) size(vol_var,2)],proj_ang2',80,[proj_ang2(1,i) proj_ang2(2,i) proj_ang2(3,i)]);
        w_proj=tom_apply_weight_function(var_stack(:,:,i),w_func);
        tom_backproj3d_euler(vol_var,w_proj,proj_ang2(1,i),proj_ang2(2,i),proj_ang2(3,i),[0 0 0]);
        if (mod(i,50)==0)
            fprintf('%s','.');
        end;
    end;
    fprintf('%s \n', ['...done! ' ]);

    all_vols(:,:,:,iiii)=vol;
    all_vars(:,:,:,iiii)=vol_var;
    all_vols_rs(:,iiii)=reshape(vol,[size(vol,1).*size(vol,2).*size(vol,3) 1]);
    
    tom_emwrite(['var_vols/rec_' num2str(iiii) '.em'],vol);
    
    clear('proj_ang2');
    
    
end;

vol_sum=sum(all_vols,4);

var_3d=tom_calc_variance(all_vols_rs');

var_3d=reshape(var_3d,[size(vol,1) size(vol,2) size(vol,3)]);


