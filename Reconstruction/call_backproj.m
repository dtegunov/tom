function call_backproj()

big_vol='big_vol.em';
tmp_folder='chunk_vols/';
proj_path='tmp_files/test*.em';    
org_size=[256 256 128]; %only tested 4 power of 2
chunk_size=[256 256 32]; %only tested 4 power of 2


warning off; mkdir(tmp_folder); warning on;
base_p=fileparts(proj_path);

if (isempty(base_p)==0)
    base_p=[base_p '/'];
end;
d=dir(proj_path);


num_of_chunks=org_size(3)./chunk_size(3);

z_off_incre=chunk_size(3);
z_off=-(z_off_incre.*(num_of_chunks./2))+(z_off_incre/2);

for ii=1:num_of_chunks
   disp(['Reconstruction chunk vol Nr: ' num2str(ii)]);
   vol_chunk=ones(chunk_size,'single');
   for i=1:length(d);
        proj_tmp=tom_emreadc([base_p d(i).name]);
        Tiltaxis=proj_tmp.Header.Tiltaxis;
        Tiltangle=proj_tmp.Header.Tiltangle;
        proj_tmp=single(proj_tmp.Value);
        tom_backproj3d(vol_chunk,proj_tmp,Tiltaxis,Tiltangle,[0 0 z_off]);
        disp(['   Backproj: ' base_p '/' d(i).name]);
    end;
    z_off=z_off+z_off_incre;
    disp(['Writing ' tmp_folder '/ch_vol_' num2str(ii) '.em']);
    tom_emwrite([tmp_folder '/ch_vol_' num2str(ii) '.em'],vol_chunk);
end;

clear('vol_chunk');
disp(['Cat chunk vols']);
vol=tom_emreadc([tmp_folder '/ch_vol_' num2str(1) '.em']);
disp(['Adding: ' tmp_folder '/ch_vol_' num2str(1) '.em']);
vol=vol.Value;
for i=2:num_of_chunks
    disp(['Adding: ' tmp_folder '/ch_vol_' num2str(i) '.em']);
    vol_tmp=tom_emreadc([tmp_folder '/ch_vol_' num2str(i) '.em']);
    vol=cat(3,vol,vol_tmp.Value);
end;
tom_emwrite(big_vol,vol);





