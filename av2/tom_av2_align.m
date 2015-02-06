function [angle_out shift_out ccc aligned_part_sum]=tom_av2_align(ref,im,mask_im,mask_cc_rot,mask_cc_trans,filter_st,num_of_iterations,demo)
%TOM_AV2_ALIGN performs iterative rotational, translational alignment of an image (im) relative to reference (ref)
%
%   [angle_out shift_out ccc aligned_part_sum]=tom_av2_align(ref,im,mask_im,mask_cc_rot,mask_cc_trans,filter_st,num_of_iterations,demo)
%
%PARAMETERS
%
%  INPUT
%   ref                 reference image
%   im                  image to be aligned
%   mask_im             mask for both, reference and image
%   mask_cc_rot         mask for polar xcorellation function, rotation alignment
%   mask_cc_trans       mask for polar xcorellation function, translation alignment
%   filter_st           filter structure:
%                        filter.Apply:  1 apply filter, 2 use default values, 0 not
%                        filter.Value:  vector with parameters [low, high, smooth ...]
%                        filter.Method: i.e. 'circ', 'quadr', 'bandpass'
%                        filter.Space: 'real' or 'fourier'
%                        filter.Times: 'apply filter
%                        n-times
%   num_of_iterations   number of iterations for iterative alignment, all
%                        rotations and shifts are added
%   demo                flag (1: on, 0: off) for demo mode via graphical interface,
%                        all alignment steps are shown
%  
%  OUTPUT
%   angle_out           rotation angle in degree, sum of all
%                        rotation angles
%   shift_out           shifts in x and y, sum of all shifts
%   ccc                 normalized xcorrelation coefficient
%   aligned_part_sum	aligned input image im
%
%EXAMPLE
%   read image
%    ref=tom_emread('/fs/pool/pool-bmsan-apps/tom_dev/data/2d/mona.em'); ref=ref.Value;
%    r=ref(40:75,35:65);
%    r=tom_norm(r,'phase');r=tom_smooth(r,8);
%    ref=tom_paste(zeros(64,64),r,[33-size(r,1)./2 33-size(r,2)./2]);
%    rotate and shift ref
%    im=tom_shift(tom_rotate(ref,45),[0 3]);
%    mask_im=tom_spheremask(ones(64),24,5);
%    mask_cc_rot=ones(128,32);
%    mask_cc_trans=ones(64,64);
%    filter.Apply=0; 
%    num_of_iterations=5;
%    demo=1;
%    [angle_out shift_out ccc aligned_part_sum]=tom_av2_align(ref,im,mask_im,mask_cc_rot,mask_cc_trans,filter,num_of_iterations,demo);
%    figure; tom_imagesc(ref);title('reference');
%    figure; tom_imagesc(im);title('image');
%    figure; tom_imagesc(aligned_part_sum);title('aligned image');
%    figure; tom_imagesc(aligned_part_sum-ref);title('difference');
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by SN/FB 01/24/06
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


%parse inputs

if (sum(size(im)==size(ref))~=2)
    error('images have different size!');
else
    im_sz=size(ref);
    middle_im=floor(im_sz./2)+1;
end;

switch nargin
    case 2
        mask_im=ones(im_sz);
        mask_cc_rot=ones(size(tom_cart2polar(ones(im_sz))));
        mask_cc_trans=ones(im_sz);
        filter_st.Apply=0;
        num_of_iterations=1;
        demo=0;
    case 3
        mask_cc_rot=ones(size(tom_cart2polar(ones(im_sz))));
        mask_cc_trans=ones(im_sz);
        filter_st.Apply=0;
        num_of_iterations=1;
        demo=0;
    case 4
        mask_cc_trans=ones(im_sz);
        filter_st.Apply=0;
        num_of_iterations=1;
        demo=0;
    case 5
        filter_st.Apply=0;
        num_of_iterations=1;
        demo=0;
    case 6
        demo=0;
        num_of_iterations=1;
    case 7
        demo=0;
    case 8

    otherwise
        error('wrong number of arguments!');
end;

if (isempty(mask_im))
    mask_im=ones(im_sz);
end;

if (isempty(mask_cc_rot))
    mask_cc_rot=ones(size(tom_cart2polar(ones(im_sz))));
end;

if (isempty(mask_cc_trans))
    mask_cc_trans=ones(im_sz);
end;

if (isempty(num_of_iterations))
    num_of_iterations=1;
end;

if (isempty(demo))
    demo=0;
end;

if (isempty(filter_st))
    filter_st.Apply=0;
end;

%demo=1;


if (demo~=0)
    demo_al(ref,'ref',0);
    demo_al(im,'im',0);
end;

im_initial=im;
ref_org=ref;

angle_vect=zeros(num_of_iterations,3);
trans_vect=zeros(num_of_iterations,3);

for i=1:num_of_iterations

    
    
    %mask and filter reference (ref) and particle (im)
    im_org=im;
    im=tom_apply_filter(im,filter_st).*mask_im;
    ref=tom_apply_filter(ref_org,filter_st).*mask_im;
    %ref=tom_apply_filter(ref_org,filter_st).*tom_spheremask(ones(160),20,5);
    
    
    
    
    % determine rotation angle by polar xcorrelation
    im_polar=tom_cart2polar(im);    
    ref_polar=tom_cart2polar(ref);
   % ccf_rot=tom_corr(im_polar,ref_polar,'norm',tom_cart2polar(mask_im));
    ccf_rot=tom_corr(im_polar,ref_polar,'norm');
    ccf_rot=ccf_rot.*mask_cc_rot;
    [ccc_pos_rot ccc]=tom_peak(ccf_rot,'spline');
    angle=(360./size(ccf_rot,2)).*((size(ccf_rot,2)./2+1)-ccc_pos_rot(2));
    im_rot=tom_rotate(im_org,-angle);
    im_rot=im_rot.*mask_im;

    % determine translation angle by cartesian xcorrelation
    %ccf_trans=tom_corr(im_rot,ref,'norm',mask_im);
    ccf_trans=tom_corr(im_rot,ref,'norm');
    ccf_trans=ccf_trans.*mask_cc_trans;
    [ccc_pos_trans ccc]=tom_peak(ccf_trans,'spline');
    shift=ccc_pos_trans-middle_im;

    if (demo~=0)
        demo_al(ref_polar,'ref_polar',0);
        demo_al(im_polar,'im_polar',0);
        demo_al(ccf_rot,'ccf_rot',ccc_pos_rot);
        demo_al(ref,'ref_rot',0);
        demo_al(im_rot,'im_rot',0);
        demo_al(ccf_trans,'ccf_trans',ccc_pos_trans);
    end;

    % store all subsequent rotations and translations
    angle_vect(i,1)=angle;
    trans_vect(i,1)=shift(1);
    trans_vect(i,2)=shift(2);
    
    % add all subsequent rotations and translations
    [angle_out shift_out]=tom_sum_rotation(-angle_vect,trans_vect,'rot_trans');
    aligned_part_sum=tom_shift(tom_rotate(im_initial,angle_out(1)),[shift_out(1) shift_out(2)]);
    
    im=aligned_part_sum;
    if (demo~=0)
        demo_al(aligned_part_sum,'im_alg',i);
    end;

    %disp('t');
    
end;

% apply full filter!!!
aligned_part_sum=tom_apply_filter(aligned_part_sum,filter_st).*mask_im;
ccf=tom_corr(aligned_part_sum,ref,'norm');
ccf=ccf.*mask_cc_trans;
[ccc_pos ccc]=tom_peak(ccf,'spline');

if (demo~=0)
    demo_al(aligned_part_sum,'im_alg_final');
end;
 
% demo options
function demo_al(im,flag,val)
if isempty(findobj('tag','demo_align2d'))
    figure;
    set(gcf,'tag','demo_align2d');
    set(gcf,'Position',[6 32 766 1092]);
    set(gcf,'Name','tom_av2_alignment_demo');
    set(gcf,'DoubleBuffer','on');
%     set(gcf,'Toolbar','none');
%     set(gcf,'MenuBar','none');
    drawnow;
else
    figure(findobj('tag','demo_align2d'));
end;

if (strcmp(flag,'im'))
    subplot(9,3,1); tom_imagesc(im,'noinfo');title('particle');axis off;
    drawnow;
end;


if (strcmp(flag,'ref'))
    subplot(9,3,3); tom_imagesc(im,'noinfo');
    title('reference');axis off;
    drawnow;
end;

if (strcmp(flag,'im_polar'))
    subplot(9,3,4); tom_imagesc(im,'noinfo');title('particle polar');axis off;
    drawnow;
end;

if (strcmp(flag,'ref_polar'))
    subplot(9,3,6); tom_imagesc(im,'noinfo');
    title('reference polar');axis off;
    drawnow;
end;

if (strcmp(flag,'ccf_rot'))
    subplot(9,3,8); tom_imagesc(im,'noinfo');title('xcf polar');axis off;
    hold on; plot(val(1),val(2),'r+'); hold off;
    drawnow;
end;

if (strcmp(flag,'im_rot'))
    subplot(9,3,10); tom_imagesc(im,'noinfo');title('particle rotated');axis off;
    drawnow;
end;

if (strcmp(flag,'ref_rot'))
    subplot(9,3,12); tom_imagesc(im,'noinfo');
    title('reference');axis off;
    drawnow;
end;

if (strcmp(flag,'ccf_trans'))
    subplot(9,3,14); tom_imagesc(im,'noinfo');title('xcf trans');axis off;
    hold on; plot(val(1),val(2),'r+'); hold off;
    drawnow;
end;

if (strcmp(flag,'im_alg'))
    if val<10
        subplot(9,3,15+val);
        tom_imagesc(im,'noinfo');
        title(['particle rot&trans, it: ' num2str(val)]);axis off;
    else
        subplot(9,3,18);
        tom_imagesc(im,'noinfo');
        title(['particle rot&trans, it: ' num2str(val)]);axis off;
    end;
    drawnow;
end;

if (strcmp(flag,'im_alg_final'))
        subplot(9,3,25);
        tom_imagesc(im,'noinfo');
        title('particle rot&trans, final');axis off;
    drawnow;
end;
if (strcmp(flag,'im_alg_sum'))
        subplot(9,3,26);
        tom_imagesc(im,'noinfo');
        title('particle rot&trans, finalsum');axis off;
    drawnow;
end;
if (strcmp(flag,'im_alg_diff'))
        subplot(9,3,27);
        tom_imagesc(im,'noinfo');
        title('particle rot&trans, difference');axis off;
    drawnow;
end;

if (strcmp(flag,'new_ref_rotated_trans'))
    subplot(9,3,21); tom_imagesc(im,'noinfo');
    title('new reference rot&trans');axis off;
    drawnow;
end;
















