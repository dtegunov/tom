function tom_av2_em_classify3d_mooped(models_path,samp_phi,samp_psi,samp_theta,num,snr,outpath,out_sel)


model_names=dir(models_path);
pathstr = fileparts(models_path);

fid=fopen(out_sel,'w');

st=tom_xmippdocread('models/model_1.doc');
zz=0;
for mod_count=1:length(model_names)
    mod=tom_emread([pathstr '/' model_names(mod_count).name]);
    mod=mod.Value;
    
    for i_phi=samp_phi(1):samp_phi(2):samp_phi(3)
        for ii_psi=samp_psi(1):samp_psi(2):samp_psi(3)
            for iii_theta=samp_theta(1):samp_theta(2):samp_theta(3)
                for i=1:num
                    zz=zz+1;
                    mod_rot=tom_rotate(mod,[i_phi ii_psi iii_theta]);
                    tmp_proj=squeeze(sum(mod_rot,2))+0.2*rand(size(mod_rot,1),size(mod_rot,2));
                    tom_spiderwrite([outpath num2str(zz) '.spi'],tmp_proj);
                    fprintf(fid,[outpath num2str(zz) '.spi 1 \n']);
                    st(zz).name=[outpath num2str(zz) '.spi'];
                    st(zz).rot=0;
                    st(zz).tilt=i_phi;
                    st(zz).psi=0;
                    st(zz).flip=0;
                    st(zz).ref=zz;
                end;
                
            end;
        end;
    end;
    
end;

st=st(1:zz);
tom_xmippdocwrite('test.doc',st);
fclose(fid);