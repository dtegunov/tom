function vol=tom_av2_apply_alignment(align2d,iter_num)
%TOM_TOM_AV2_APPLY_ALIGNMENT creates ...
%
%   vol=tom_av2_apply_alignment(align2d,iter_num)
%
%PARAMETERS
%
%  INPUT
%   align2d             ...
%   iter_num          	...
%  
%  OUTPUT
%   vol                 ...
%
%EXAMPLE
%   ... = tom_av2_apply_alignment(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by FB mm/dd/yy
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

 vol=back_proj(align2d,iter_num);




function vol=back_proj(align2d,iter_num)

%transfer often used data
angular_scan_eu=align2d(iter_num,1).angular_scan_euler;
angular_scan=align2d(iter_num,1).angular_scan;
size_stack=align2d(iter_num,1).stack_size;
size_vol=align2d(iter_num,1).model.size;
size_part=align2d(iter_num,1).model.size_part;
shift_corr_flag=align2d(iter_num,1).shiftcorr;

vol=zeros(size_vol,'single'); %allocate some memory 


size_stack(3)=650;


%find used angleclasses ...build theta vector for weighting 
%future rip ...weightig according 2 num of paricles per class could be performed here !
zz=1;
for i=1:size_stack(3)
%for i=1:36
     eu_4_weight(zz,:)=align2d(iter_num,i).angleclass.angle_euler;
     zz=zz+1;
end;
% eu_4_weight(:,1)=0;
% eu_4_weight(:,2)=0;


%theta_4_weight=angular_scan_eu(:,3);

zz=1;
%write out weighting function ...speed up!
% for i=1:size(angular_scan_eu,1)
%     eu=angular_scan_eu(i,:);
%     if (align2d(iter_num,i).avg.num_array(i)~=0)
%         if (size_part(1)~=size_part(2))
%             thickn=size_part(1).*abs(cos((angular_scan(1,i).*(pi./180)))) + size_part(2).*abs(sin((angular_scan(1,i).*(pi./180)))); %generic sample thickness for barrels
%         else
%             thickn=size_part(1);
%         end;
%         w=tom_weight2d_exact([size_stack(1) size_stack(2)],eu,eu_4_weight,thickn);
%         tom_emwrite(['tmp/weighting_' num2str(i) '.em'],w);
%         zz=zz+1;
%         disp(zz);
%     end;
% end;


mask=tom_spheremask(ones(size_stack(1)),size_stack(1)./2-6,3,[size_stack(1)./2+1 size_stack(1)./2+1 1]);
%[w_func_a,c_mask]= create_masks(size(vol,1)./2-4,[size(vol,1) size(vol,1)]);

% perform weighting and backprojection
%for i=1:36
 for i=1:size_stack(3)
    eu=align2d(iter_num,i).angleclass.angle_euler;
    shiftxy=align2d(iter_num,i).angleclass.shiftxy;
    
    %read particle
    part=tom_emreadc(align2d(iter_num,1).filename,'subregion',[0 0 i],[size_stack(1)-1 size_stack(2)-1 0]);
    part=double(part.Value);

    
    
%     [a,b,iso]=tom_calc_pixel_thresh(-tom_norm(part,1),5000,0.01);
%     
%     mask=zeros(size(part));
%     mask=tom_paste(mask,ones(120,50),[20 55]);
%     mask=tom_filter(mask,6);
%     part=mask.*-(iso~=0);
    
    %correct xy-shift
    if (shift_corr_flag==1)
        part=tom_shift(part,shiftxy);
    end;
    
    %load weighting function and apply it 
   % w_func=tom_emreadc(['tmp/weighting_' num2str(align2d(iter_num,i).angleclass.proj_nr) '.em']);
%     [w_func]= create_masks([160 160]);
%     w_func=tom_emheader(w_func);
    
   % part=tom_rotate(part,90,'linear');
  
    
    %part = real(ifft2((fft2(part).*fftshift(w_func.Value) ))); 
%     if (isempty(find(isnan(w_func.Value)))==0)
%          disp('nan w_func');
%     end;
    %part=tom_rotate(part,-90,'linear');
    if (isempty(find(isnan(part)))==0)
        disp('nan part');
    end;
    
    %part=tom_weight3d_euler(part,size(vol,1)./2-2,eu(1),eu(2),theta_4_weight',i,thickn);
    vol=single(vol);
    num_of_p= align2d(iter_num,1).avg.num_array(align2d(iter_num,i).angleclass.proj_nr);
    part=part./num_of_p;
    part=single(part.*mask);
    tom_backproj3d_euler(vol,part,eu(1),eu(2),eu(3),[0 0 0]);
   
    
    
    if (isempty(find(isnan(vol)))==0)
        disp('nan vol');
    end;
    
    disp(i);
    %command window print
%     store.i=i;
%     store=disp_out(store,'estimate_time');
%     store=disp_out(store,'progress');
end;
eu_4_weight=align2d(1,1).angular_scan_euler
[w]=tom_weight3d_exact([160 160 160],eu_4_weight,160);
vol = real(ifftn((fftn(vol).*fftshift(w) ))); 


 if (isempty(find(isnan(vol)))==0)
        disp('removing nan');
%         vol=tom_remove_nan(vol);
 end;

 %apply symmetry
 
% vol=tom_symref(vol,2,'C2',[180 0 161]);
 
 vol=vol.*tom_spheremask(ones(size(vol)),size(vol,1)./2-6,3);
%  [a b vol]=tom_calc_isosurface(-vol,2400,3.4,0.01);
%  xxx= ones(size(vol));
%  yyy = tom_cylindermask(xxx,24,2,(size(vol)./2)+1);
%  vol=tom_rotate(yyy,[270 90 90]).*vol~=0.*tom_spheremask(ones(size(vol)),64,2);
%  vol=-vol;
 function [w_func,c_mask]= create_masks(dim_h)

[xh,yh] = ndgrid(-dim_h(1)./2:((dim_h(2)./2)-1));
w_func = tom_norm(abs(yh),1);
