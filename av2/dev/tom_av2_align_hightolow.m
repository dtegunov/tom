function tom_av2_align_hightolow(inalignfile,outalignfile,Filter,Mask_trans,Mask_rot,verboseflag)
% TOM_AV2_ALIGN_HIGHTOLOW automatically picks particles from the low
% defocus series by using the particles from the high defocus series
%
% Syntax:
% tom_av2_align_hightolow(inalignfile,outalignfile,Filter,Mask_trans,Mask_rot,verboseflag)
%
% Input:
% inalignfile:   filename of the input alignment file
% outalignfile:  filename of the output alignment file
% Filter:  
% Mask_trans:   
% Mask_rot:      
% verboseflag:   1 for output messages, 0 for quiet operation (optional)
%
% created: 19/1/06 AK
%
%   Copyright (c) 2006
%   TOM toolbox for Electron Tomography
%   Max-Planck-Institute for Biochemistry
%   Dept. Molecular Structural Biology
%   82152 Martinsried, Germany
%   http://www.biochem.mpg.de/tom


try
    s = load(inalignfile);
catch
    error('Input alignment file could not be loaded.');
end


inalign = s.align2d;

filename_cache = '';

outalign = inalign;
h = tom_reademheader(inalign(1).filename);
h.Header.Size = h.Header.Size./4;


if (nargin==2)
    Filter.Apply=2;
    Mask_trans.Apply=2;
    Mask_trans.Value=[h.Header.Size(1) h.Header.Size(2)];
    Mask_rot.Apply=0;
    Mask_rot.Value=[h.Header.Size(1)./2 h.Header.Size(2).*2];
    verboseflag = 0;
end;

if (nargin==3)
    Filter.Apply=2;
    Mask_trans.Apply=2;
    Mask_trans.Value=[h.Header.Size(1) h.Header.Size(2)];
    Mask_rot.Apply=0;
    Mask_rot.Value=[h.Header.Size(1) h.Header.Size(2)];
    verboseflag = 0;
end;

if (nargin==4)
    Mask_rot.Apply=0;
    Mask_rot.Value=[h.Header.Size(1) h.Header.Size(2)];
    verboseflag = 0;
end;

if (nargin==5)
    verboseflag = 0;
end;

if (isempty(Filter)==1)
     Filter.Apply=2;
end;


if (isempty(Mask_rot)==1)
    Mask_rot.Apply=0;
    Mask_rot.Value=[h.Header.Size(1)./2 h.Header.Size(2).*2];
end;

if (isempty(Mask_trans)==1)
    Mask_trans.Apply=0;
    Mask_trans.Value=[h.Header.Size(1) h.Header.Size(2)];
end;

mask_im = tom_spheremask(ones(h.Header.Size(1),h.Header.Size(2)),h.Header.Size(2)./2-20,16);
mask_rot=tom_create_mask(Mask_rot);
mask_trans=tom_create_mask(Mask_trans);

for i=1:size(inalign,2)
    outalign(i).filename = regexprep(inalign(i).filename,'high', 'low');
    
    %only align files if particle is on a new file
    if ~isequal(filename_cache,inalign(i).filename)
        im_h = tom_emreadc(inalign(i).filename,'resample',[4 4 1]);
        im_l = tom_emreadc(outalign(i).filename,'resample',[4 4 1]);
        im_h = tom_norm(double(im_h.Value),'phase');
        im_l = tom_norm(double(im_l.Value),'phase');
    
        [angle_out shift_out aligned_part]=tom_av2_align(im_h,im_l,mask_im,mask_rot,mask_trans,Filter,2,0);
        shift_out=-shift_out;
        r = tom_pointrotate([shift_out(1).*4, shift_out(2).*4,0],angle_out(1),0,0);
        filename_cache = inalign(i).filename;
    end

    outalign(i).position.x = round(inalign(i).position.x + r(1));
    outalign(i).position.y = round(inalign(i).position.y + r(2));

    if verboseflag == 1 & mod(i,10) == 0 
        disp([num2str(i) ' of ' num2str(size(inalign,2)) ' particles done']);
    end
    
end

align2d = outalign;
save(outalignfile,'align2d');

if verboseflag == 1
    disp('Finished');
end