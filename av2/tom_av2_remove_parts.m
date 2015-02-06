function img=tom_av2_remove_parts(align2d,flag,reduce,radius,bin,im)
%TOM_AV2_REMOVE_PARTS remove particlels in the picklist form the image
%
%   img=tom_av2_remove_parts(align2d,flag,reduce,im)
%
%PARAMETERS
%
%  INPUT
%   align2d             align2d struct
%   flag                how the part coord should be replaced  
%                       (localNoise) localMean or perm
%   reduce              reduce value for the replcement patch in %./100 
%   radius              replacement radius refering 2 unbinned image
%   bin                 (0) binning of the coordinates if no image given also
%                        of the image
%  
%   im                  (opt.) image where the particles should be removed 
%                       ...use for speed
%
%  OUTPUT
%    img                image with replaced parts  
%
%
%EXAMPLE
%
%  new=tom_av2_remove_parts(pl_red,'localNoise',-0.01,76,2,imb);
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by SN/FB 01/24/06
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

if (nargin<5)
    im='';
end;
    

if (isempty(im))
    im=tom_emread(align2d(1,1).filename);
    im=tom_bin(im.Value,bin);
end;

img=im;

rad=radius./(2.^bin);

mask=tom_spheremask(ones(rad.*2,rad.*2),rad);
idx=find(mask==1);

img=im;

for i=1:size(align2d,2)
    coord=round([align2d(1,i).position.x align2d(1,i).position.y]./(2.^bin));
    org=tom_cut_out(im,[coord-rad],[rad.*2 rad.*2]);
    mea=mean(org(:))-(reduce.*mean(org(:)));
    st_d=std(org(:));
    new=org;
    
    if (strcmp(flag,'localNoise'))
        noise=(tom_norm(rand(size(new)),'mean0+1std').*st_d)+mea;    
        new=(new.*(mask==0))+(noise.*mask);
    end;
    if (strcmp(flag,'localMean'))
        noise=(tom_norm(ones(size(new)),'mean0+1std').*st_d)+mea;    
        new=(new.*(mask==0))+(noise.*mask);
    end;
    
    img=tom_paste(img,new,[coord-rad]);
    
  
end;












