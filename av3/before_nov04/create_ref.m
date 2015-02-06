pdbdata = tom_pdbread('GroEL.pdb');
emmap = tom_pdb2em(pdbdata,4.1 , 128);
cm = tom_cm(emmap);
cent = [size(emmap,1)/2+1 size(emmap,2)/2+1 size(emmap,3)/2+1];
cm = cm - cent; clear cent;
emmap = tom_shift(emmap,-cm);clear cm;%center
% if necessary ...
%[toi,eigvec,eigval] = tom_moi(emmap)
ctf = -tom_create_ctf(-12,emmap, 0.41, 300, 2, 0.3);
%design low-pass
[x,y,z]=meshgrid( -floor(size(ctf,1)/2):floor(size(ctf,1)/2)-1,...
    -floor(size(ctf,2)/2):floor(size(ctf,2)/2)-1, -floor(size(ctf,3)/2):floor(size(ctf,3)/2)-1);
lowp = sqrt(x.^2 +y.^2+z.^2) <= 15; 
clear y; clear z;clear x;
% filter after 2nd zero 
ctf_filt= ctf.*lowp; clear lowp;
femmap = ifftn(ifftshift(fftshift(fftn(emmap)).*ctf_filt));
%adjust pixel-size: 0.265 -> 1.7nm
%fribo_ref = fribo_ref(101-15:101+15,101-15:101+15,101-15:101+15);
%ribo_ref = real(ifftn(ifftshift(fribo_ref)));
femmap = tom_bin(femmap,1);
tom_dspcub(femmap);
%confine to 21 pixels
%ribo_ref = ribo_ref(6:26,6:26,6:26);
tom_emwrite('GroEL.rem',femmap);
femmap = tom_emheader(femmap);
