function tom_fit_ctf_local(align2d,box_size,diff2global,corr_flag,mtf,cut_off,output_align,output_stack,output_param)

%
% tom_fit_ctf_local(align2d,1024,2,'flip',[],400,'/fs/sun11/lv01/pool/pool-nickell2/Titan/040310/04032010/log/xmipp_rec/4Nickell2/align2d_local.mat', ...
% '/fs/sun11/lv01/pool/pool-nickell2/Titan/040310/04032010/log/xmipp_rec/4Nickell2/stack_local.em','em_stack')
% tom_fit_ctf_local(align2d,1024,1,'flip',[],400,'/fs/sun11/lv01/pool/pool-nickell2/Titan/040310/04032010/log/xmipp_rec/4Nickell2/align2d_local_1mue.mat', ...
% '/fs/sun11/lv01/pool/pool-nickell2/Titan/040310/04032010/log/xmipp_rec/4Nickell2/stack_local_1mue.em','em_stack')


if (exist('output_param','var')==0)
    output_param='em_stack';
end;

if (exist('mtf','var')==0)
    mtf='';
end;



if (strcmp(output_param,'em_stack'))
    stack_corr=zeros(align2d(1,1).radius.*2,align2d(1,1).radius.*2,size(align2d,2));
    align2d_out=align2d;
    stack_count=1;
end;


%get unique list of filenames
for i=1:size(align2d,2)
    all_names{i}=align2d(1,i).filename;
end;
name_list=unique(all_names);

h=tom_reademheader(name_list{1});
sz=h.Header.Size;

diff2global=diff2global*1e-06;



%parfor i=1:length(name_list)
 
fp=fopen([output_align '.txt'],'w');
fclose(fp);

parfor i=1:length(name_list)

    result{i}=fit_and_correct_local(name_list{i},align2d,all_names,diff2global,box_size,sz,corr_flag,cut_off,mtf,[output_align '.txt']);
end;

all=1;
for i=1:length(name_list)
    all=cat(1,all,result{i}.done);
end;
all=all(2:end);


stack=zeros(align2d(1,1).radius.*2,align2d(1,1).radius.*2,length(find(all==1)));

zz=1;
for i=1:length(name_list)
    for ii=1:length(result{i}.done)
        if (result{i}.done(ii)==1) 
            stack(:,:,zz)=result{i}.stack_corr(:,:,ii);
            align2d(1,zz)=result{i}.align2d(1,ii);
            zz=zz+1;
        end;
    end;
end;

%cut alignment file
tmp=align2d(1,1:zz-1);
clear('align2d');
align2d=tmp;

save(output_align,'align2d');
tom_emwrite(output_stack,stack);




function result=fit_and_correct_local(name,align2d,all_names,diff2global,box_size,sz,corr_flag,cut_off,mtf,log_name)

   
   fp=fopen(log_name,'a');
 
    [a b c]=fileparts(name);
    idx=find(ismember(all_names,name));
    
    stack_corr=zeros(align2d(1,1).radius.*2,align2d(1,1).radius.*2,length(idx));
    align2d_out=align2d(1,length(idx));
    done=zeros(length(idx),1);
    
    try
        st_out=my_load([name '.mat']);
    catch
        disp('Error no global fit!!');
        fprintf(fp,['Error: Cannot open '  name '.mat \n']);
    end;
    
    %check search parameters
    dz_d=(st_out.Search.Dz_search(2)-st_out.Search.Dz_search(1))./10;
    dz=st_out.Fit.Dz_det;
    st_out.Search.Dz_search=(dz-(diff2global/2)):dz_d:(dz+(diff2global/2));
    st_out.Search.Dz_delta_search=st_out.Fit.Dz_det;
    st_out.Search.Phi_0_search=st_out.Fit.Phi_0_det;
    fclose(fp);
    
    for ii=1:length(idx)
        try
            [pos shift]=get_position([align2d(1,idx(ii)).position.x align2d(1,idx(ii)).position.y],[box_size box_size],[sz(1) sz(2)]);
            img_box = tom_emreadc3(name, [pos(1),pos(2),0 box_size,box_size,1],[1,1,1], [1,1,1]);
            img_box=img_box.Value;
            
            fit_out=fit_image(st_out,img_box);
            disp([name ' Dz global: ' num2str(dz) ' local: ' num2str(fit_out.Dz_det) ' diff: ' num2str(fit_out.Dz_det-dz) ]);
            fp=fopen(log_name,'a');
            fprintf(fp,[name ' ; ' num2str(dz) ' ; ' num2str(fit_out.Dz_det) ' ; ' num2str(fit_out.Dz_det-dz) ' \n']);
            fclose(fp);
            img_box=corr_image(img_box,fit_out,corr_flag,cut_off,mtf);
            
            stack_corr(:,:,ii)=tom_cut_out(tom_move(img_box,-shift),'center',[align2d(1,1).radius.*2 align2d(1,1).radius.*2]);
            align2d_out(1,ii)=align2d(1,idx(ii));
            % align2d_out(1,ii).local_def=fit_out.Dz_det;
            done(ii)=1;
        catch
            disp('error in cutting particle! ...continue');
            fp=fopen(log_name,'a');
            fprintf(fp,['Error: Cannot process '  name ' \n']);
            fclose(fp);
        end;
     end;

%fclose(fp);
     
result.stack_corr=stack_corr;
result.align2d=align2d_out;
result.done=done;


function st_out=my_load(name)

load(name);



function [pos_out shift]=get_position(pos_center,box_size,img_size)

box_size=ceil(box_size./2);

%calc box pos
x_min=pos_center(1)-box_size(1);
y_min=pos_center(2)-box_size(2);

x_max=pos_center(1)+box_size(1);
y_max=pos_center(2)+box_size(2);


pos_out=pos_center-box_size-1;

pos_real=pos_out;

%check bounderies

if (x_min<1) 
    pos_out(1)=1;
end;

if (y_min<1)
    pos_out(2)=1;
end;    

if (x_max > img_size(1)) 
    pos_out(1)=pos_out(1) - (x_max-img_size(1))-1;
end;

if (y_max > img_size(2))
    pos_out(2)=pos_out(2) - (y_max-img_size(2))-1;
end;    

shift=pos_real-pos_out;


function Fit=fit_image(ctf_struct,in)


img_size= ctf_struct.Search.ps_size;
mask_in_radius=ctf_struct.Fit.mask_inner_radius;
mask_out_radius=ctf_struct.Fit.mask_outer_radius;
mask_in = tom_spheremask(ones(img_size),mask_in_radius,0,[img_size(1)./2+1 img_size(2)./2+1 1]);
mask_out = tom_spheremask(ones(img_size),mask_out_radius,0,[img_size(1)./2+1 img_size(2)./2+1 1]);
mask=mask_out-mask_in;


EM=ctf_struct.Fit.EM;


inner_rad_bg=ctf_struct.Search.mask_inner_radius_bg;
outer_rad_bg=ctf_struct.Search.mask_outer_radius_bg;
Search_tmp=ctf_struct.Search;



ps=tom_calc_periodogram(double(in),img_size(1),img_size(1)/16);
ps=(log(fftshift(ps)));

[decay decay_image]=calc_decay(ps,inner_rad_bg,outer_rad_bg,32);
background_corrected_ps=double(ps-decay_image);

ps=background_corrected_ps.*mask;
warning off;
[Fit]=tom_fit_ctf(ps,EM,Search_tmp);
warning on;



function image=corr_image(image,fit,method,corr_cut_off,mtf)

    
image=tom_correct_for_ctf_and_mtf_new(double(image),fit,method,corr_cut_off,mtf);




function pos_out=write_image(image,name,flag)




% stack_corr=stack_corr(:,:,find(done==1));
% align2d_out=align2d_out(1,find(done==1));







%for i=1:length(name_list)
% for i=1:10
% diff2global
%     [a b c]=fileparts(name_list{i});
%     idx=find(ismember(all_names,name_list{i}));
%     img_name=name_list{i};
%     st_out=my_load([name_list{i} '.mat']);
%     %check search parameters
%     dz_d=st_out.Search.Dz_search(2)-st_out.Search.Dz_search(1);
%     dz=st_out.Fit.Dz_det;
%     st_out.Search.Dz_search=(dz-(diff2global./2)):dz_d:(dz+(diff2global./2));
%     st_out.Search.Dz_delta_search=st_out.Fit.Dz_det;
%     st_out.Search.Phi_0_search=st_out.Fit.Phi_0_det;
%     
%     for ii=1:length(idx)
%         [pos shift]=get_position([align2d(1,idx(ii)).position.x align2d(1,idx(ii)).position.y],[box_size box_size],[sz(1) sz(2)]);
%         img_box = tom_emreadc3(img_name, [pos(1),pos(2),0 box_size,box_size,1],[1,1,1], [1,1,1]);
%         img_box=img_box.Value;
%         
%         fit_out=fit_image(st_out,img_box);
%         img_box=corr_image(img_box,fit_out,corr_flag,cut_off,mtf);
%         if (strcmp(output_param,'em_stack'))
%             stack_corr(:,:,i)=tom_cut_out(tom_move(img_box,-shift),'center',[align2d(1,1).radius.*2 align2d(1,1).radius.*2]);
%             align2d_out(1,idx(ii))=align2d(1,idx(ii));
%             align2d_out(1,i).local_def=fit_out.Dz_det;
%             done(idx(ii))=1;
%         else
%             %pos_out=write_image(image,name,flag);
%         end;
%     end;
%     
%     disp(name_list{i});
%     
% end;
% 
% stack_corr=stack_corr(:,:,find(done==1));
% align2d_out=align2d_out(1,find(done==1));
% 
% disp(' ');
% 















