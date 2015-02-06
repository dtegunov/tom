function [vol_var,vol_avg,vol_var_even,vol_var_odd]=tom_xmippvariancemap3d_penczek(doc_filename,num_of_models,proj_per_model,filter_kernel,mask,tmp_filepath,outputlevel)
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
%    [vol_var,vol_avg,vol_var_even,vol_var_odd]=tom_xmippvariancemap3d_penczek('26S_out__it00012doc');
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
    tmp_filepath=[pwd '/rec_files4var'];
end;

if nargin < 7
    outputlevel=1;
end;


if exist(tmp_filepath,'dir')==0
    mkdir(tmp_filepath);
    unix(['chmod ugo+rwx ' tmp_filepath]);
end;

tic;

%vol_var=zeros(size(im_tmp.Value,1),size(im_tmp.Value,1),size(im_tmp.Value,1));

if (outputlevel==1); fprintf('%s ', 'Building Volumes: ' ); end;

num_of_pix=size(im_tmp.Value,1).^3;
parfor i=1:num_of_models
    
    model=rec_model(st_all,proj_per_model,outputlevel);
    %norm volumes 2 mean 1 due debub and numerics for calc of variance
    
    if (isempty(find(isnan(model))  )==0 )
       disp([ num2str((length(find(isnan(model)))./num_of_pix).*100) ' % nan ' num2str(length(find(isnan(model)))) ' found !' ]);
       model=tom_remove_nan(model,mean(mean2(model(find(isnan(model)==0)) )  ) ); 
    end;
    model=tom_norm(model,'mean0+1std')+1;
    model=tom_xraycorrect3d_local(double(model));
    
    tom_emwrite([tmp_filepath '/rec_' num2str(i) '.em'],model);
    
    if (outputlevel==1); fprintf('%s','.'); end;
    
 end;

toc;


[vol_var vol_avg vol_var_even vol_var_odd]=tom_av3_calc_variance(tmp_filepath,num_of_models,mask,filter_kernel,outputlevel);
 
 
  
 
    function vol=rec_model(st_all,proj_per_model,outputlevel)
        
        
        
        for iiii=1:10
            
            
            
            if (outputlevel==2); disp(['Building ' num2str(iiii) ' Volume']); end;
            
            idx=round(length(st_all).*rand(proj_per_model,1));
            disp('check out for random doublicate bug');
            idx=idx+(idx==0);
            st=st_all(idx);
            
            tmp=zeros(size(st,2),1);
            proj_angles=zeros(size(st,2),3);
            tmp_angles=zeros(size(st,2),3);
            
            
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
            
            proj_ang2=zeros(3,num_of_ref);
            
            im_tmp=tom_spiderread([ st(1).name]);
            stack=zeros(size(im_tmp.Value,1),size(im_tmp.Value,2),num_of_ref);
            stack2=zeros(size(im_tmp.Value,1),size(im_tmp.Value,2),num_of_ref);
            count=zeros(num_of_ref);
            count2=zeros(num_of_ref);
            class=zeros(size(st,2),1);
           
            if (outputlevel==2); fprintf('%s \n', ['...done! ' ]); end;
            
            if (outputlevel==2); fprintf('%s ', ['Building Classes: ' ]); end;
            
            for i=1:size(st,2)
                
                im=tom_spiderread([ st(i).name]);
                im.Value=tom_xraycorrect2(im.Value);
                if (st(i).flip==0)
                    im_tmp_alg=tom_rotate(tom_move(im.Value,round([st(i).xoff st(i).yoff])),tmp_angles(i,3));
                    stack(:,:,st(i).ref)=stack(:,:,st(i).ref)+im_tmp_alg;
                    count(st(i).ref)=count(st(i).ref)+1;
                    [aa tttemp]=tom_eulerconvert_xmipp(st(i).rot, st(i).tilt, 0);
                    proj_ang2(:,st(i).ref)=tttemp;
                    class(i)=st(i).ref;
                else
                    im.Value=tom_mirror(im.Value,'x');
                    im.Value=im.Value;
                    im_tmp_alg=tom_rotate(tom_move(im.Value,round([-st(i).xoff st(i).yoff])),tmp_angles(i,3));
                    stack2(:,:,st(i).ref)=stack2(:,:,st(i).ref)+im_tmp_alg;
                    count2(st(i).ref)=count2(st(i).ref)+1;
                    [aa tttemp]=tom_eulerconvert_xmipp(st(i).rot, st(i).tilt, 0);
                    proj_ang2(:,st(i).ref)=tttemp;
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
           
             stack=stack+stack2;
            
            zz=1;
            for i=1:size(stack,3)
                if (std2(stack(:,:,i))~=0)
                    proj_ang_clean(:,zz)=proj_ang2(:,i); 
                    zz=zz+1;
                end;
            end;
           
            
            mask_sp=tom_spheremask(ones(size(stack,1),size(stack,1),size(stack,1)),round(size(stack,1)./2)-4,3);
            %mask_sp=tom_spheremask(ones(size(stack,1),size(stack,1)),round(size(stack,1)./2)-10,3 );
           
            if (outputlevel==2); fprintf('%s ', ['Reconstructing: ' ]); end;
            %backproject image stack
           
            vol=zeros(size(stack,1),size(stack,1),size(stack,1),'single');
            for i=1:size(stack,3)
                if (std2(stack(:,:,i))~=0)
                   % w_func=tom_calc_weight_functionc([size(vol,1) size(vol,2)],proj_ang_clean',size(vol,1),[proj_ang2(1,i) proj_ang2(2,i) proj_ang2(3,i)]);
                   % w_proj=tom_apply_weight_function(stack(:,:,i),w_func,mask_sp);
                    %danger Hack fb
                    w_proj=stack(:,:,i);
                    
                    tom_backproj3d_euler(vol,w_proj,proj_ang2(1,i),proj_ang2(2,i),proj_ang2(3,i),[0 0 0]);
                    if (mod(i,50)==0)
                        if (outputlevel==2); fprintf('%s','.'); end;
                    end;
                end;
            end;
            %weight it in 3d baby!
%             
            w_func=tom_calc_weight_functionc(size(vol),proj_ang_clean',size(vol,1));
            
            inter=w_func(33,32,:);
            inter=inter+w_func(32,33,:);
            inter=inter+w_func(34,33,:);
            inter=inter+w_func(33,34,:);
            inter=inter./4;
            
            w_func(33,33,:)=inter;
            vol=tom_apply_weight_function(vol,w_func,mask_sp);
            
            
            if (outputlevel==2); fprintf('%s \n', ['...done! ' ]); end;
            
             if (isempty(find(isnan(vol))   ))
                break;
             end;
            
        end;
        
        
        
        
        
        
        
        
function I=tom_xraycorrect3d_local(varargin)


if nargin<1
    error('Not enough Input Arguments');
    return;
end
if nargin==1
    J=varargin{1};
    I=J;
    s=2;
elseif nargin==2
    J=varargin{1};
    I=J;
    s=varargin{2};
else
    error('Too many input Arguments');
    return;
end

d = std2(J(:));
meanimg = mean(J(:));
[x] = find(J>(meanimg+18.*d) | J<(meanimg-18.*d)); % factor 10 added by SN
[x,y,z] = ind2sub(size(J), x);
if isempty(x)
    return;
else
    for i=1:size(x,1)
        if ((x(i)-s > 0) && (y(i)-s > 0) && (z(i)-s > 0) && ...
            (x(i)+s <= size(I,1)) && (y(i)+s <= size(I,2)) && (z(i)+s <= size(I,3)) )
            I(x(i),y(i),z(i))=mask(I,x(i),y(i),z(i),d,meanimg,s);
        else
            I(x(i),y(i),z(i))=mean(I(:));
        end;
    end
end
if size(x,1)>0
    disp(['values outside range found: ' num2str(size(x,1))]); %changed by SN
end;


%****** SubFunction mask ***********

function c=mask(A,xcoor,ycoor,zcoor,dev,meanimg,s)

[xdim ydim zdim]=size(A);
a=A(xcoor-s:xcoor+s,ycoor-s:ycoor+s,zcoor-s:zcoor+s);
t=find(a>(meanimg+dev) | a<(meanimg-dev));
if isempty(t)
    c=sum(sum(a))/((2*s+1)*(2*s+1));
elseif (size(a,1)*size(a,2)*size(a,3) == size(t,1))%bug fixed FF
    c=0;
else
    a(t)=0;
    c=sum(a(:))/(((2*s+1)*(2*s+1)*(2*s+1))-size(t,1));
end


 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
  
 
 
 

