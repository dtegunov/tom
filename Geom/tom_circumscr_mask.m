function mask=tom_circumscr_mask(bin_vol,shape,offset)
%  TOM_CIRCUMSCR_MASK circumscribes mask with a given 
%  
%     vol=tom_spheremask(vol, radius,sigma,center);
%  
%  PARAMETERS
%  
%    INPUT
%     bin_vol             binary vol
%     shape               shape of the mask (cyl,sphere ..)
%     offset              offest 2 boarders
%    
%    OUTPUT
%     vol                 mask
%  
%  EXAMPLE
%     
%  
%  REFERENCES
%  
%  SEE ALSO
%     ...
%  
%     
%  
%     Nickell et al., 'TOM software toolbox: acquisition and analysis for electron tomography',
%     Journal of Structural Biology, 149 (2005), 227-234.
%  
%     Copyright (c) 2004-2007
%     TOM toolbox for Electron Tomography
%     Max-Planck-Institute of Biochemistry
%     Dept. Molecular Structural Biology
%     82152 Martinsried, Germany
%     http://www.biochem.mpg.de/tom

sz=size(bin_vol);

if (strcmp(shape,'cyl'))
    sig=find(squeeze(sum(sum((bin_vol),2),1))>0);
    mi=min(sig)-offset(1)-1; ma=max(sig)+offset(1)+1;  
    sig=find(squeeze(sum(sum((bin_vol),3),1))>0);
    mit=min(sig); mat=max(sig); rad=mat-mit; 
    mask=tom_cylindermask(ones(sz),round(rad./2)+offset(2));
    mask(:,:,1:mi)=0; mask(:,:,ma:end)=0;
end;





