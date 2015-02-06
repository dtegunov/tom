function [Fit]=tom_fit_ctf(in,EM,Search)

%
% in=tom_emreadc(['carbon_1.em']);
% ps=tom_calc_periodogram_parallel(single(in.Value),256,0,256./16);

% tom_show_ctf_fit(psb3_clean_background, Dz_det, Dz_delta_det, Phi_0_det, decay_det,3.637e-10,3e5,2e-3,15,64,0);
% 

Dz_search=Search.Dz_search;
Dz_delta_search=Search.Dz_delta_search;
Phi_0_search=Search.Phi_0_search;
decay_search=Search.decay_search;
amplitude_contrast_search=Search.amplitude_contrast_search;
mask_in_radius=Search.mask_inner_radius;
mask_out_radius=Search.mask_outer_radius;
decay_part_coh_ill_search=Search.decay_part_coh_ill_search;
decay_energy_spread_search=Search.decay_energy_spread_search;

pix_size=EM.Objectpixelsize;
voltage=EM.Voltage;
Cs=EM.Cs;
Cc=EM.Cc;



img_size=size(in);

% create the 2d function
Ny = 1./(2.0*pix_size);
[x y] = meshgrid(-Ny:2*Ny/img_size(1):Ny-Ny/img_size(1), -Ny:2*Ny/img_size(2):Ny-Ny/img_size(2));
k = sqrt(x.^2.0+y.^2.0);
k_2 = x.^2.0+y.^2.0;

% create the masks and the mask indices.
mask_in = tom_spheremask(ones(img_size),mask_in_radius,0,[img_size(1)./2+1 img_size(2)./2+1 1]);
mask_out = tom_spheremask(ones(img_size),mask_out_radius,0,[img_size(1)./2+1 img_size(2)./2+1 1]);
mask=mask_out-mask_in;
ind_mask=find(mask);
ind_mask_size=size(ind_mask,1);

 %in_norm=mask.*tom_norm(in,'mean0+1std',mask);
in_norm_index=(in(ind_mask)-mean(in(ind_mask)))./std(in(ind_mask));

% pre-allocate the output variables
corr_all=zeros(size(Dz_search,2),size(decay_search,2),size(Dz_delta_search,2),size(Phi_0_search,2),size(amplitude_contrast_search,2),size(decay_part_coh_ill_search,2),size(decay_energy_spread_search,2));

% pre-calculate the constants
voltage_relativistic=voltage*(1+voltage/1022000.0); % relativistic
lambda=sqrt(150.4/voltage_relativistic)*10^-10; % wavelength relativistic
Phi_0_search=Phi_0_search./180.*pi; % convert degree in rad
pi_2=pi./2;
term_1=lambda.^3.*Cs.*k.^4;
term_1_ind=term_1(ind_mask);
at2=atan2(x,y);
at2_ind=at2(ind_mask);
lambda_2=2.*lambda;
Dz_delta_search_2=Dz_delta_search./2;
decay_search_Ny=decay_search.*Ny;
k_ind=k(ind_mask);
k_2_ind=k_2(ind_mask);
Cc_voltage=Cc./voltage;

% for later - improved parallelisation: idx=1;for i=1:size(Dz_search,2);for ii=1:size(Dz_delta_search,2);Dz_Dz_delta_search(idx,1)=Dz_search(i);Dz_Dz_delta_search(idx,2)=Dz_delta_search_2(ii);idx=idx+1;end;end;

% loop over all four variables
%parfor Dz_run=1:size(Dz_search,2)
for Dz_run=1:size(Dz_search,2)
    Dz_test=Dz_search(Dz_run);
    corr_tmp=inner_loop(Dz_test,Dz_delta_search_2,Phi_0_search,decay_search_Ny,pi_2,term_1_ind,lambda_2,at2_ind,k_ind,k_2_ind,in_norm_index,amplitude_contrast_search,decay_part_coh_ill_search,Cs,Cc_voltage,decay_energy_spread_search);
    corr_all(Dz_run,:,:,:,:,:,:,:)=corr_tmp;
end;

% apply normalization;
% find maximum and corresponding indices
[i1 i2 i3 i4 i5 i6 i7]=ind2sub(size(corr_all),find(corr_all==max(max(max(max(max(max(max(corr_all))))))),1));

Fit.Dz_det=Dz_search(i1(1));
Fit.Dz_delta_det=Dz_delta_search(i3(1));
Fit.Phi_0_det=Phi_0_search(i4(1))./pi.*180;
Fit.amplitude_contrast_det=amplitude_contrast_search(i5(1));
Fit.decay_det=decay_search(i2(1));
Fit.decay_part_coh_ill_det=decay_part_coh_ill_search(i6(1));
Fit.decay_energy_spread_det=decay_energy_spread_search(i7(1));
Fit.corr_all=corr_all./(ind_mask_size-1);
Fit.mask_inner_radius=Search.mask_inner_radius;
Fit.mask_outer_radius=Search.mask_outer_radius;
Fit.EM=EM;


function corr_tmp=inner_loop(Dz_test,Dz_delta_search_2,Phi_0_search,decay_search_Ny,pi_2,term_1_ind,lambda_2,at2_ind,k_ind,k_2_ind,in_norm_index,amplitude_contrast_search,decay_part_coh_ill_search,Cs,Cc_voltage,decay_energy_spread_search)

% pre-allocate the output variables
corr_tmp=zeros(size(decay_search_Ny,2),size(Dz_delta_search_2,2),size(Phi_0_search,2));

for Dz_delta_run=1:size(Dz_delta_search_2,2)
    Dz_delta_test=Dz_delta_search_2(Dz_delta_run);
    for decay_run=1:size(decay_search_Ny,2)
        decay_test=decay_search_Ny(decay_run);
        if decay_test == 0
            decay_test=ones(size(k_ind));
        else
            decay_test = exp(-(k_ind./(decay_test)).^2);
        end;
        for Phi_0_run=1:size(Phi_0_search,2)
            Phi_0_test=Phi_0_search(Phi_0_run);
            for amplitude_contrast_run=1:size(amplitude_contrast_search,2)
                amplitude_contrast_test=amplitude_contrast_search(amplitude_contrast_run);
                for decay_part_coh_ill_run=1:size(decay_part_coh_ill_search,2)
                    decay_part_coh_ill_test=decay_part_coh_ill_search(decay_part_coh_ill_run);
                    
                    if decay_part_coh_ill_test == 0
                        E_i_k=1;
                    else
                        q_0=decay_part_coh_ill_test./((Cs.*(lambda_2./2).^3).^0.25); % see formula 2.19a
                        E_i_k=exp((-pi.^2.*q_0.^2.*(Cs.*(lambda_2./2).^3.*k_ind.^3-Dz_test.*lambda_2./2.*k_ind).^2)); % for a gaussian source distribution
%                        decay_part_coh_ill_test=E_i_k;
                    end;
                    
                    for decay_energy_spread_run=1:size(decay_energy_spread_search,2)
                        decay_energy_spread_test=decay_energy_spread_search(decay_energy_spread_run);
                        
                        if decay_energy_spread_test == 0
                            E_e_k=1;
                        else                            
                            delta_z = Cc_voltage.*decay_energy_spread_test;
                            E_e_k=exp(-(pi.*delta_z.*lambda_2.*2.*(k_ind.^2)./2)); % for a gaussian source distribution                            
%                            decay_energy_spread_test=E_e_k;
                        end;
                        
                        w=pi_2.*(term_1_ind-lambda_2.*(Dz_test+Dz_delta_test.*sin(2.0.*(at2_ind-Phi_0_test))).*k_2_ind);
                        amplitude=cos(w);
                        phase    =sin(w);
                        % better readable
                        % amplitude=cos(pi/2.0.*(lambda.^3.*Cs.*k.^4-2.*lambda.*(Dz_test+(Dz_delta_test./2.0).*sin(2.0.*(atan2(x,y)-Phi_0_test))).*k.^2));
                        % phase=sin(pi/2.0.*(lambda.^3.*Cs.*k.^4-2.*lambda.*(Dz_test+(Dz_delta_test./2.0).*sin(2.0.*(atan2(x,y)-Phi_0_test))).*k.^2));
                        %              amplitude_norm=-mask.*decay.*tom_nor
                        %              m((amplitude.*amplitude),'mean0+1std',mask);
                        amplitude_index_2=E_e_k.*E_i_k.*decay_test.* ...
                            ((sqrt(1-amplitude_contrast_test.^2).*phase+amplitude_contrast_test.*amplitude).^2);
                        amplitude_index_2=(amplitude_index_2-mean(amplitude_index_2))./std(amplitude_index_2);
                        corr_tmp(decay_run,Dz_delta_run,Phi_0_run,amplitude_contrast_run,decay_part_coh_ill_run,decay_energy_spread_run)=sum(amplitude_index_2.*in_norm_index);
                    end;                    
                end;                
            end;
        end;
    end;
end;
