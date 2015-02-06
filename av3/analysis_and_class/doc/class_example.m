%th_open = tom_pdbread('../../../templates/x-ray/thermosome/hsp_mod_a6d_open.pdb');
%th_closed = tom_pdbread('../../../templates/x-ray/thermosome/hsp_model_closed.pdb');
%th_open_em = tom_pdb2em(th_open, 4, 64);
%th_closed_em = tom_pdb2em(th_closed, 4, 64);
%th_open_fil = tom_filter(th_open_em,1);
%th_closed_fil = tom_filter(th_closed_em,1);
%tom_emwrite('th_closed.em',th_closed_fil);
%tom_emwrite('th_open.em',th_open_fil);
%th_open_fil = tom_emread('th_open.em');th_open_fil = th_open_fil.Value;
% added for correct sampling
%th_open_fil = tom_bandpass(th_open_fil,0,14,3);%th_open_fil = tom_filter(th_open_fil,1.5,'circ');
%th_closed_fil = tom_emread('th_closed.em');th_closed_fil = th_closed_fil.Value;
%th_closed_fil = tom_bandpass(th_closed_fil,0,14,3);%th_closed_fil = tom_filter(th_closed_fil,1.5,'circ');
%ctf = -tom_create_ctf(-3, th_open_fil, 0.4, 300, 2, 0.4);
tic
mask = tom_spheremask(ones(64,64,64),29,1);
%[mv mn mx stv sigma] = tom_dev(th_open_fil);
for miswedge=0:10:40
    for iSNR=0:2
        figure(1);clf;
        SNR = 1/(10^iSNR);
        motl = tom_emread('../motl_test.em');motl = motl.Value;
        wedgelist(1,1) = 1;
        wedgelist(2,1) = -(90-miswedge);wedgelist(3,1) = 90-miswedge;
        name = ['../test_' 'SNR' num2str(SNR) 'wedge' num2str(miswedge)];
        [ccc conorm] = tom_pca(motl, name, wedgelist,mask,0,15,1,0);
        name = ['./ccc_' 'SNR' num2str(SNR) 'wedge' num2str(miswedge) '.em'];
        tom_emwrite(name,ccc);
        name = ['./conorm_' 'SNR' num2str(SNR) 'wedge' num2str(miswedge) '.em'];
        tom_emwrite(name,conorm);
        [coeff, latent, explained] = av3_pcacov(ccc);
        %[vec val] = eigs(ccc,6);
        %for ieig=1:6
        %    ncoeff(:,ieig) = sqrt(latent(ieig))*coeff(:,ieig);
        %end;
        % now incorporated into av3_pcacov
        ncoeff = coeff(:,1:6);%FF 05/13/05
        idx = av3_kmeans(ncoeff(:,2:3), 2);
        % design pca-stack
        pcastack = ncoeff;
        pcastack(:,7) = idx;% last column: classes
        pcastack(21,1:6) = latent(1:6);pcastack(21,7) = trace(ccc);
        pcastack(22,1:6) = explained(1:6);
        name = ['./cpcastack_' 'SNR' num2str(SNR) 'wedge' num2str(miswedge) '.em'];
        tom_emwrite(name,pcastack);
        idx1 = find(idx == 1);idx2 = find(idx == 1);
        figure(2);clf;
        plot(ncoeff(idx1,2),ncoeff(idx1,3),'+');hold on;
        plot(sqrt(latent(2))*coeff(idx2,2),sqrt(latent(3))*coeff(idx2,3),'o');
        name = ['./coeff_' 'SNR' num2str(SNR) 'wedge' num2str(miswedge) '.em'];
        tom_emwrite(name,coeff(:,1:6));
        name = ['./eigvals_' 'SNR' num2str(SNR) 'wedge' num2str(miswedge) '.em'];
        tom_emwrite(name,latent(1:6));
        name = ['./eigvals_rel_' 'SNR' num2str(SNR) 'wedge' num2str(miswedge) '.em'];
        tom_emwrite(name,explained(1:6));
        %[eigstack, evstack] = av3_prep_class_in_EM(ccc,6,name,motl,wedgelist,0,8);
        %name = ['eigstack_' 'SNR' num2str(SNR) 'wedge' num2str(miswedge) '.em'];
        %tom_emwrite(name,eigstack);
        %name = ['evstack_' 'SNR' num2str(SNR) 'wedge' num2str(miswedge) '.em'];
        %tom_emwrite(name,evstack);
        %proj_c= zeros(size(motl,2),6);
        for ieig=1:6
            name = ['../test_' 'SNR' num2str(SNR) 'wedge' num2str(miswedge)];
            average = av3_eigvec2vol(motl, ncoeff(:,ieig), name,wedgelist,0);%rec eigenvector
            name = ['./eigvec_' num2str(ieig) 'SNR' num2str(SNR) 'wedge' num2str(miswedge) '.em'];
            tom_emwrite(name, average);
            %proj_c(:,ieig) = av3_constr_ccc(average, motl, './test_noise', wedgelist,mask,0,20,0,0);
        end;        
        %tom_emwrite(['reconst_SNR' num2str(SNR) 'wedge' num2str(miswedge) '.em'],val(ieig,ieig)*proj_c);
        %tom_emwrite(['test_SNR' num2str(SNR) 'wedge' num2str(miswedge) '_stack.em'],noisystackforem);
    end;
end;
toc
%Elapsed time is 13854.441026 seconds.