function [im_handle h_amplitude h_phase defocus_label_points_x defocus_label_points_y amplitude_norm]=tom_show_ctf_fit(in, EM, Fit)

% tom_show_ctf_fit(tom_filter(wow,3), Dz_det, Dz_delta_det, Phi_0_det, decay_det,3.6e-10,3e5,2e-3,8,96,0.07)


pix_size=EM.Objectpixelsize;
voltage=EM.Voltage;
Cs=EM.Cs;
Cc=EM.Cc;

Dz_det=Fit.Dz_det;
Dz_delta_det=Fit.Dz_delta_det;
Phi_0_det=Fit.Phi_0_det;
amplitude_contrast=Fit.amplitude_contrast_det;
decay_det=Fit.decay_det;
decay_part_coh_ill_det=Fit.decay_part_coh_ill_det;
decay_energy_spread_det=Fit.decay_energy_spread_det;
%Fit.corr_all=corr_all./(ind_mask_size-1);
mask_in_radius=Fit.mask_inner_radius;
mask_out_radius=Fit.mask_outer_radius;

img_size=size(in);
% create the 2d function
Ny = 1./(2.0*pix_size);
[x y] = meshgrid(-Ny:2*Ny/img_size(1):Ny-Ny/img_size(1), -Ny:2*Ny/img_size(2):Ny-Ny/img_size(2));
k = sqrt(x.^2.0+y.^2.0);

% create the masks and the mask indices.
mask_in = tom_spheremask(ones(img_size),mask_in_radius,0,[img_size(1)./2+1 img_size(2)./2+1 1]);
mask_out = tom_spheremask(ones(img_size),mask_out_radius,0,[img_size(1)./2+1 img_size(2)./2+1 1]);
mask=mask_out-mask_in;

in_norm=mask.*tom_norm(in,'mean0+1std',mask);

% pre-calculate the constants
voltage_relativistic=voltage*(1+voltage/1022000.0); % relativistic
lambda=sqrt(150.4/voltage_relativistic)*10^-10; % wavelength relativistic
Phi_0_det=Phi_0_det./180.*pi; % convert degree in rad
pi_2=pi./2;
term_1=lambda.^3.*Cs.*k.^4;
at2=atan2(x,y);
lambda_2=2.*lambda;




decay_test=decay_det;
if decay_test == 0
    decay=1;
else
    decay = exp(-(k./(decay_test.*Ny)).^2);
end;
Dz_test=Dz_det;
Dz_delta_test=Dz_delta_det;
Phi_0_test=Phi_0_det;
amplitude=cos(pi_2.*(term_1-lambda_2.*(Dz_test+(Dz_delta_test./2.0).*sin(2.0.*(at2-Phi_0_test))).*k.^2));
phase=sin(pi_2.*(term_1-lambda_2.*(Dz_test+(Dz_delta_test./2.0).*sin(2.0.*(at2-Phi_0_test))).*k.^2));
% better readable
%amplitude=cos(pi/2.0.*(lambda.^3.*Cs.*k.^4-2.*lambda.*(Dz_test+(Dz_delta_test./2.0).*sin(2.0.*(atan2(x,y)-Phi_0_test))).*k.^2));
%phase=sin(pi/2.0.*(lambda.^3.*Cs.*k.^4-2.*lambda.*(Dz_test+(Dz_delta_test./2.0).*sin(2.0.*(atan2(x,y)-Phi_0_test))).*k.^2));

if decay_part_coh_ill_det == 0
    decay_part_coh_ill_det=ones(size(k));
    E_i_k=1;
else
    q_0=decay_part_coh_ill_det./((Cs.*(lambda_2./2).^3).^0.25); % see formula 2.19a
    E_i_k=exp((-pi.^2.*q_0.^2.*(Cs.*(lambda_2./2).^3.*k.^3-Dz_test.*lambda_2./2.*k).^2)); % for a gaussian source distribution
end;
if decay_energy_spread_det == 0
    decay_energy_spread_det=ones(size(k));
    E_e_k=1;
else
    delta_z = Cc.*decay_energy_spread_det./voltage;
    E_e_k=exp(-(pi.*delta_z.*lambda.*k.^2./2)); % for a gaussian source distribution
end;

amplitude_norm=mask.*tom_norm(E_e_k.*E_i_k.*decay.*((sqrt(1-amplitude_contrast.^2).*phase + amplitude_contrast.*amplitude).^2),'mean0+1std',mask);
merged=zeros(img_size);
merged(img_size(1)./2+1:end,:)=amplitude_norm(img_size(1)./2+1:end,:);
merged(1:img_size(1)./2+1,:)=in_norm(1:img_size(1)./2+1,:);
merged_norm=mask.*tom_norm(merged,'mean0+1std',mask);

rot_phase=tom_rotate(amplitude_norm',Fit.Phi_0_det-45).*mask;
%rot_phase=tom_rotate(phase',Fit.Phi_0_det-45).*mask;
[c h_phase_rot]=contour(rot_phase,1);
[x1 x2]=find(c(2,:)==img_size(1)./2+1);  v=c(1,x2); v=v(find((v-(img_size(1)./2+1+mask_in_radius+1))>0)); 
v=v(find(v<((img_size(1)./2+1+mask_out_radius-1)))); v2=ones(size(v)).*img_size(1)./2+1; % plot(v,v2,'.');
[y1 y2]=find(c(1,:)==img_size(1)./2+1);  vy=c(2,y2); vy=vy(find((vy-(img_size(1)./2+1+mask_in_radius+1))>0));
vy=vy(find(vy<((img_size(1)./2+1+mask_out_radius-1)))); vy2=ones(size(v)).*img_size(1)./2+1; % plot(vy2,vy,'r.');

for i=1:size(v,2); 
    try
        tmp=tom_pointrotate2d([v(i) v2(i)],(Fit.Phi_0_det-45),size(phase)./2);
        vxr(i,1)=tmp(1);
        vxr(i,2)=tmp(2);
    catch
    end;
end;

for i=1:size(vy,2); 
    try
        tmp=tom_pointrotate2d([vy2(i) vy(i)],(Fit.Phi_0_det-45),size(phase)./2); 
        vyr(i,1)=tmp(1); 
        vyr(i,2)=tmp(2);
    catch
    end;
end;

defocus_label_points_x=vxr;
defocus_label_points_y=vyr;


im_handle=imagesc(merged_norm');axis image; colormap gray;
hold on;
[c h_amplitude]=contour(amplitude',1);
set(h_amplitude,'color','green');

%[c h_phase]=contour(phase',1);
[c h_phase]=contour((sqrt(1-amplitude_contrast.^2).*phase + amplitude_contrast.*amplitude)',1);
set(h_phase,'color','red');
%[x1 x2]=find(c(2,:)==129);  v=c(1,x2); v2=c(2,x2); plot(v,v2,'.');
hold off;
drawnow;



%title(['Defocus: ' num2str(Dz_det) ' Dz-delta: ' num2str(Dz_delta_det) ' Phi-0: ' num2str(Phi_0_det.*180./pi) ' Decay: ' num2str(decay_det)]);
%figure; plot(merged_norm(:,129))

%mytexstr = '$\cos(\frac{\Pi}{2} (\lambda^3 C_s k^4 - 2 \lambda (Dz+\frac{\Delta Dz}{2} \sin ( 2 ( \arctan(x,y) - \phi_0 ))) k^2 ))$';
%h = text('string',mytexstr,'interpreter','latex','fontsize',35,'units','norm','pos',[.05 .5]);
