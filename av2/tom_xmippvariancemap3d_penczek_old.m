function [vol_var,vol_avg,vol_var_even,vol_var_odd]=tom_xmippvariancemap3d_penczek(doc_filename,num_of_models,proj_per_model,filter_kernel,mask,tmp_filepath,ouputlevel)
%tom_xmippvariancemap3d_penczek calculates a 3d variance Map 
%
%   [vol_var,var_stack,vol,stack]=tom_xmippvariancemap3d_penczek(doc_filename)
%
%   calculates 3d variance volume by the Penzcek bootstrapping method.
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
%    [vol_var,var_stack,vol,stack]=tom_xmippvariancemap3d_penczek('26S_out__it00012doc');
%
%REFERENCES
%
%   Penczek Pawel A; Yang Chao; Frank Joachim; Spahn Christian M T
%   Estimation of variance in single-particle reconstruction using the bootstrap technique.
%   Journal of structural biology 2006;154(2):168-83.
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
im_tmp=tom_spiderread(st_all(1).name);


if nargin < 2
    num_of_models=500;
end;

if nargin < 3
    proj_per_model=5000;
end;

if nargin < 4
    filter_kernel=1;
end;

if nargin < 5
    mask=ones(size(im_tmp.Value,1),size(im_tmp.Value,1),size(im_tmp.Value,1));
end;

if nargin < 6
    tmp_filepath='rec_files4var';
end;

if nargin < 7
    outputlevel=1;
end;


if exist(tmp_filepath,'dir')==0
    mkdir(tmp_filepath);
end;


vol_avg=zeros(size(im_tmp.Value,1),size(im_tmp.Value,1),size(im_tmp.Value,1));
vol_var=zeros(size(im_tmp.Value,1),size(im_tmp.Value,1),size(im_tmp.Value,1));

if (outputlevel==1); fprintf('%s ', 'Building Volumes: ' ); end;

for iiii=1:num_of_models

    if (outputlevel==2); disp(['Building ' num2str(iiii) ' Volume']); end;
    
    idx=round(length(st_all).*rand(proj_per_model,1));
    idx=idx+(idx==0);
    st=st_all(idx);


   if (outputlevel==2); fprintf('%s ', ['Converting Angles: ' ]); end;
    for i=1:size(st,2)
        tmp(i)=st(i).ref;
        [xx,angles] = tom_eulerconvert_xmipp(st(i).rot, st(i).tilt, st(i).psi);
        proj_angles(i,:) = angles';
        tmp_angles(i,:)=[st(i).rot st(i).tilt st(i).psi];
        if (mod(i,2000)==0)
           if (outputlevel==2); fprintf('%s','.'); end;
        end;
    end;
    if (outputlevel==2); fprintf('%s \n','done!'); end;


    if (outputlevel==2); fprintf('%s ', ['Allocating Memory: ' ]); end;
    num_of_ref=max(tmp);
    %im_tmp=tom_spiderread(['./' st(1).name]);
    im_tmp=tom_spiderread([ st(1).name]);
    stack=zeros(size(im_tmp.Value,1),size(im_tmp.Value,2),num_of_ref);
    stack2=zeros(size(im_tmp.Value,1),size(im_tmp.Value,2),num_of_ref);
    count=zeros(num_of_ref);
    count2=zeros(num_of_ref);
    class=zeros(size(st,2),1);
    feat_stack=zeros(size(st,2),size(im_tmp.Value,1).*size(im_tmp.Value,1));
    if (outputlevel==2); fprintf('%s \n', ['...done! ' ]); end;

    if (outputlevel==2); fprintf('%s ', ['Building Classes: ' ]); end;
    for i=1:size(st,2)
        
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
           if (outputlevel==2); fprintf('%s','.'); end;
        end;
    end;
    if (outputlevel==2); fprintf('%s \n', ['...done! ' ]); end;

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

    if (outputlevel==2); fprintf('%s ', ['Reconstructing: ' ]); end;
    %backproject image stack
    stack=stack+stack2;
    vol=zeros(size(stack,1),size(stack,1),size(stack,1),'single');
    for i=1:size(proj_ang2,2)
        w_func=tom_calc_weight_functionc([size(vol,1) size(vol,2)],proj_ang2',80,[proj_ang2(1,i) proj_ang2(2,i) proj_ang2(3,i)]);
        w_proj=tom_apply_weight_function(stack(:,:,i),w_func);
        tom_backproj3d_euler(vol,w_proj,proj_ang2(1,i),proj_ang2(2,i),proj_ang2(3,i),[0 0 0]);
        if (mod(i,50)==0)
            if (outputlevel==2); fprintf('%s','.'); end;
        end;
    end;
    if (outputlevel==2); fprintf('%s \n', ['...done! ' ]); end;
 

    tom_emwrite([tmp_filepath '/rec_' num2str(iiii) '.em'],vol);
    
    clear('proj_ang2');
    
    vol_avg=vol_avg+tom_norm(double(vol),'mean0+1std');
   
    if (outputlevel==1); fprintf('%s','.'); end;
    
 end;

 if (outputlevel==1); fprintf('%s \n','done!'); end;
 
 vol_avg=vol_avg./iiii;
 


  fprintf('%s ', ['Calculating Variance: ' ]);
  im_big=zeros((size(im_tmp.Value,1).*size(im_tmp.Value,1).*size(im_tmp.Value,1)),iiii);
  
  for i=1:iiii
     tmp=tom_emread([tmp_filepath '/rec_' num2str(i) '.em']);
     tmp=double(tmp.Value);
     tmp=tom_filter(tmp,filter_kernel);
     tmp=tmp.*mask;
     tmp=tom_norm(tmp,'mean0+1std',mask);
     %vol_var=vol_var+((tmp-vol_avg).*(tmp-vol_avg));
     im_big(:,i)=reshape(tmp,[size(im_tmp.Value,1).*size(im_tmp.Value,1).*size(im_tmp.Value,1) 1]);
     fprintf('%s','.');
  end;
  
  fprintf('%s \n', ['done! ' ]);
  
 pause(1);
 % vol_var=vol_var./iiii;
 vol_var=tom_calc_variance(im_big');
 vol_var=reshape(vol_var,[size(im_tmp.Value,1) size(im_tmp.Value,1) size(im_tmp.Value,1)]);
  
 
 idx=1:iiii;
 idx_even=find(mod(idx,2).*idx);
 idx_odd=find((mod(idx,2)==0).*idx);
 
 vol_var_even=tom_calc_variance(im_big(:,idx_even)');
 vol_var_even=reshape(vol_var_even,[size(im_tmp.Value,1) size(im_tmp.Value,1) size(im_tmp.Value,1)]);
  
 
 vol_var_odd=tom_calc_variance(im_big(:,idx_odd)');
 vol_var_odd=reshape(vol_var_odd,[size(im_tmp.Value,1) size(im_tmp.Value,1) size(im_tmp.Value,1)]);
 
 pause(1);
 
 cc=tom_corr(vol_var_even,vol_var_odd,'norm');
 [a b]=tom_peak(cc);
 
 disp(['Odd/Even Variance correlation: ' num2str(b)]);
 
 
 
 
 
 
  
  
 
 
 

