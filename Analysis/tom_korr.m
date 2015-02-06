function [corr_pic_norm]=tom_korr(pic_A,pic_B,corr_flag);

%this function computes the cross-correlation function for 2 pictures
% 
% INPUT:    picA: picture1
%           picB: picture2
%           corr_flag: determinse whether xcf,mcf or pcf is used 
%  OUTPUT:  cross-correlation function 
%  
%
%  22/11/03 SN, 24/11/03 tested and bug fixed FB
%
%   Copyright (c) 2004
%   TOM toolbox for Electron Tomography
%   Max-Planck-Institute for Biochemistry
%   Dept. Molecular Structural Biology
%   82152 Martinsried, Germany
%   http://www.biochem.mpg.de/tom



%norm the Values
meanA=mean(mean(pic_A));
meanB=mean(mean(pic_B));
ttA=pic_A-meanA;
ttB=pic_B-meanB;

% check Hack
%pic_A=ttA;
% end check Hack

ttA=ttA.*ttA;
ttB=ttB.*ttB;
ttA=sqrt(sum(sum(ttA)));
ttB=sqrt(sum(sum(ttB)));

%transform the Data
fft_dataA=(fft2(pic_A));
fft_dataB=(fft2(pic_B));


% kill middle Value
%fft_dataA(1,1)=0;
%fft_dataB(1,1)=0;

%make cross correlation
if (strcmp(corr_flag,'xcf'))
    corr_fft=(conj(fft_dataA).*(fft_dataB));
elseif (strcmp(corr_flag,'mcf'))
    corr_fft=( conj(fft_dataA).*(fft_dataB) )./(sqrt((abs(conj(fft_dataA).*(fft_dataB) ))));
elseif (strcmp(corr_flag,'pcf'))
    corr_fft=(conj(fft_dataA).*(fft_dataB))./(abs(conj(fft_dataA).*(fft_dataB)));
end;

corr_pic=real(fftshift(ifft2(corr_fft)));

%norm the correlation funktion
corr_pic_norm=corr_pic/(ttA*ttB);
















