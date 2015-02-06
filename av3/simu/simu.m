%ribo = tom_emread('Y:\templates\Ribosomes\volFEG.em');
% correct scaling: 2.93 A
%lowp = round(2.93/15*size(ribo.Value,1)); % filter to 1.5 nm res.
%ribo = tom_bandpass(ribo.Value,0,lowp,3);
%ribo = tom_bin(ribo);
%tribo = tom_filter(ribo.Value,8);
%[mn mxv minv] = tom_dev(tribo);ribo = ribo.Value-mn;
%ribo = ribo/(mxv-mn)*0.6;ribo = ribo+1;
% assume SNR proj =0.1 => SNR tomo app 1/20 *snr_proj
snr1=1;snr2=10;
dthe = 10;
proj = tom_spheremask(sum(ribo,3),30,1);
[mn mxv minv stv pvari] = tom_dev(proj);
sribo = tom_spheremask(ribo,30,1);
wedge=tom_wedge(ribo,30);
wribo = tom_spheremask(tom_fourier(ifftshift(fftshift(tom_fourier(ribo)).*wedge)),30,1);
[mn mxv minv stv vari] = tom_dev(ribo);
for ithe=0:dthe:180
    ribo_rot = tom_rotate(ribo,[0 0 ithe]);
    cc3(int8(ithe/dthe + 1)) = tom_ccc(sribo,tom_spheremask(ribo_rot,30,1),'norm');
    cc3_err1(int8(ithe/dthe + 1)) = tom_ccc(sribo,tom_spheremask(tom_error(ribo_rot,'G',0,snr1*vari),30,1),'norm');%
    cc3_err2(int8(ithe/dthe + 1)) = tom_ccc(sribo,tom_spheremask(tom_error(ribo_rot,'G',0,snr2*vari),30,1),'norm');
    cc2(int8(ithe/dthe + 1)) = tom_ccc(tom_spheremask(sum(ribo_rot,3),30,1), proj,'norm');
    cc2_err1(int8(ithe/dthe + 1)) = tom_ccc(tom_spheremask(tom_error(sum(ribo_rot,3),'G',0,pvari*snr1),30,1), proj,'norm');
    cc2_err2(int8(ithe/dthe + 1)) = tom_ccc(tom_spheremask(tom_error(sum(ribo_rot,3),'G',0,pvari*snr2),30,1), proj,'norm');
    wribo_rot = tom_spheremask(tom_fourier(ifftshift(fftshift(tom_fourier(ribo_rot)).*wedge)),30,1);
    wribo_rot_err1 = tom_spheremask(tom_fourier(ifftshift(fftshift(tom_fourier(tom_error(ribo_rot,'G',0,vari*snr1))).*wedge)),30,1);
    wribo_rot_err2 = tom_spheremask(tom_fourier(ifftshift(fftshift(tom_fourier(tom_error(ribo_rot,'G',0,vari*snr2))).*wedge)),30,1);
    cc3w(int8(ithe/dthe + 1)) = tom_ccc(wribo,wribo_rot,'norm');
    cc3w_err1(int8(ithe/dthe + 1)) = tom_ccc(wribo,wribo_rot_err1,'norm');
    cc3w_err2(int8(ithe/dthe + 1)) = tom_ccc(wribo,wribo_rot_err2,'norm');
end;
plot(0:dthe:180,cc3);hold on;plot(0:dthe:180,cc3w,'k');plot(0:dthe:180,cc2,'r');
plot(0:dthe:180,cc3_err1,':');hold on;plot(0:dthe:180,cc3w_err1,'k:');plot(0:dthe:180,cc2_err1,'r:');
plot(0:dthe:180,cc3_err2,'--');hold on;plot(0:dthe:180,cc3w_err2,'k--');plot(0:dthe:180,cc2_err2,'r--');