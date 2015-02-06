function xmipp_wien_scrt()

%variables will be adapted by awk
in_doc='xxx_doc_xxx';
df_goup_nr=xxx_999_xxx;


all_defocus=[];



%pix_size=Fit.EM.Objectpixelsize;
pix_size=pix_size.*1e-10;
voltage=Fit.EM.Voltage;
%img_size=sz_im;
sigma=Fit.decay_det;
q_0_r=Fit.decay_part_coh_ill_det;
Cs=Fit.EM.Cs;
Cc=Fit.EM.Cc;
deltaE=Fit.decay_energy_spread_det;
amplitude_contrast=Fit.amplitude_contrast_det;


ctf=tom_ctf;

ctf_theory=ctf3d(ctf,cutoff,ctf);


wiener_f = ctf_theory ./ ((ctf_theory.*ctf_theory) + wiener_const);





function w_func=ctf3d(ctf,cutoff,sx)


for ii=1:2.*sx
    for jj=1:sx
        decay_restore_3d_sphere(:,ii,jj) = decay_restore;
    end;
end;
decay_restore_3d = tom_sph2cart(decay_restore_3d_sphere);

mask=tom_spheremask(ones(sx,sx,sx),cutoff(1),0);

w_func=decay_restore_3d.*mask;