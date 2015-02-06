function [parts_corrected D Q]=phase_flip_090622_all_19595amn_allgood(Fit_start,Search)
parts=tom_emreadc3('090622_all_19595amn_allgood.em');
parts_corrected=parts;
load 090622_all_19595amn_allgood.mat

    img_size=[256 256];
    mask_in_radius=Fit_start.mask_inner_radius;
    mask_out_radius=Fit_start.mask_outer_radius;
    mask_in = tom_spheremask(ones(img_size),mask_in_radius,0,[img_size(1)./2+1 img_size(2)./2+1 1]);
    mask_out = tom_spheremask(ones(img_size),mask_out_radius,0,[img_size(1)./2+1 img_size(2)./2+1 1]);
    mask=mask_out-mask_in;

D=0;
Q=0;
i=1;
ii=0;
filename_old='';
EM=Fit_start.EM;

mtf=load('/fs/sandy01/lv03/pool/bmsan/apps/tom_dev/data/mtfs/old_mtfs/mtf_Eagle_SN109_200kV.mat');
mtf=mtf.mtf_Eagle_SN109_200kV;

tic;
%for i=1:1:2000
for i=1:1:3

    
    
    if ~isequal(filename_old,align2d(i).filename)
    in=tom_emreadc3(align2d(i).filename);

    ps=tom_calc_periodogram(double(in.Value),[256]);
    ps=(log(fftshift(ps)));
    
    [decay decay_image]=calc_decay(ps,Fit_start.decay_mask_inner_radius,Fit_start.decay_mask_outer_radius,32);
    background_corrected_ps=double(ps-decay_image);

    ps=background_corrected_ps.*mask;
    [Fit]=tom_fit_ctf(ps,EM,Search);
    end;
    ii=ii+1;    
    D(ii)=Fit.Dz_det;
    [v w]=tom_peak(Fit.corr_all);
    Q(ii)=w;
    
    [corrected_particle]=tom_correct_for_ctf_and_mtf_new(double(parts.Value(:,:,i)),Fit,'flip&mtf',40,mtf);
%    [corrected_particle]=tom_correct_for_ctf_and_mtf_new(double(parts.Value(:,:,i)),Fit,'flip',128,mtf);

    parts_corrected.Value(:,:,i)=corrected_particle;
    filename_old=align2d(i).filename;
    if ~mod(i,50)
        ii
        toc
        tic;
    end;
    
end;


