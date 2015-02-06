function tom_av2_align_hightolow(inalignfile,outalignfile,Filter,Mask_trans,Mask_rot,verboseflag)
%TOM_AV2_ALIGN_HIGHTOLOW automatically picks particles from the low ...
%
%   tom_av2_align_hightolow(inalignfile,outalignfile,Filter,Mask_trans,Mask_rot,verboseflag)
%
%  TOM_AV2_ALIGN_HIGHTOLOW automatically picks particles from the low
%  defocus series by using the particles from the high defocus series
%
%PARAMETERS
%
%  INPUT
%   inalignfile         filename of the input alignment file
%   outalignfile        filename of the output alignment file
%   Filter              ...
%   Mask_trans          ...
%   Mask_rot            ...
%   verboseflag         1 for output messages, 0 for quiet operation (optional)
%  
%  OUTPUT
%
%EXAMPLE
%   tom_av2_align_hightolow(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by AK 1/19/06
%
%   Nickell et al., 'TOM software toolbox: acquisition and analysis for electron tomography',
%   Journal of Structural Biology, 149 (2005), 227-234.
%
%   Copyright (c) 2004-2007
%   TOM toolbox for Electron Tomography
%   Max-Planck-Institute of Biochemistry
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
h.Header.Size = h.Header.Size./8;
h.Header.Size =[1024 1024];
sz = h.Header.Size(1);
center=floor(sz./2)+1;

if (nargin==2)
    Filter.Apply=1;
    Filter.Type = 'kernel';
    Filter.Value.radius=3;
    Filter.Value.method='circ';
    Filter.Value.space='real';
    Filter.Value.times=1;
    Mask_trans.Apply=1;
    Mask_trans.Type = 'sphere2d';
    Mask_trans.Value.size_x = h.Header.Size(1);
    Mask_trans.Value.size_y = h.Header.Size(2);
    Mask_trans.Value.radius=50;
    Mask_trans.Value.sigma=10;
    Mask_trans.Value.center_x=round(min(sz)./2)+1;
    Mask_trans.Value.center_y=round(min(sz)./2)+1;
    %    Mask_trans.Value=[h.Header.Size(1) h.Header.Size(2)];
    Mask_rot.Apply=1;
    Mask_rot.Type = 'sphere2d';
    Mask_rot.Value.size_x = h.Header.Size(1)./2;
    Mask_rot.Value.size_y = h.Header.Size(2).*2;
    Mask_rot.Value.radius=h.Header.Size(1)./2;
    Mask_rot.Value.sigma=10;
    Mask_rot.Value.center_x=Mask_rot.Value.size_x./2+1;
    Mask_rot.Value.center_y=Mask_rot.Value.size_y./2+1;
    verboseflag = 1;
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

%for i=1:100
for i=1:size(inalign,2)

    %inalign(i).filename = regexprep(inalign(i).filename,'high_bin1', 'high');
    outalign(i).filename = regexprep(inalign(i).filename,'high', 'low');
    
    %only align files if particle is on a new file
    if ~isequal(filename_cache,inalign(i).filename)
        im_h = tom_emreadc(inalign(i).filename,'binning',[3 3 1]);
        im_h.Value = tom_xraycorrect(im_h.Value);
        im_l = tom_emreadc(outalign(i).filename,'binning',[3 3 1]);
        im_l.Value = tom_xraycorrect(im_l.Value);
        im_h = tom_norm(double(im_h.Value),'phase');
        im_l = tom_norm(double(im_l.Value),'phase');
    
        [angle_out shift_out]=tom_av2_align(im_h,im_l,mask_im,mask_rot,mask_trans,Filter,3,0);
        filename_cache = inalign(i).filename;
        disp(inalign(i).filename);
    end

        pos_tmp=[inalign(i).position.x inalign(i).position.y]; 
        pos_tmp=pos_tmp-shift_out(1:2)'.*8;
        pos_tmp = tom_pointrotate2d(pos_tmp,-angle_out(1),[4097 4097]);
    
        outalign(i).position.x = round(pos_tmp(1));
        outalign(i).position.y = round(pos_tmp(2));

    if verboseflag == 1 && mod(i,10) == 0 
        disp([num2str(i) ' of ' num2str(size(inalign,2)) ' particles done']);
    end
   
end

align2d = outalign;
save(outalignfile,'align2d');

if verboseflag == 1
    disp('Finished');
end