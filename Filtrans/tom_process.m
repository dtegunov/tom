function out=tom_process(in,struct)
%TOM_PROCESS creates ...
%
%   out=tom_process(in,struct)
%
%PARAMETERS
%
%  INPUT
%   in                  ...
%   struct              ...
%  
%  OUTPUT
%   out                 ...
%
%EXAMPLE
%   ... = tom_process(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by ... (author date)
%   updated by ...
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


progr=struct.program;

if (struct.stack==1)
    stop=size(in,3);
end;

for ii=1:stop

if (struct.stack==1)
    im=in(:,:,ii);  
else
    im=in;
end;
    
for i=1:1000

   [field_name progr]=strtok(progr);
   
   type=getfield(getfield(struct,field_name),'Type');
   
   switch type
       case {'Filter'}
         im=tom_apply_filter(im,getfield(struct,field_name));
       case {'Mask'}
        im=im.*tom_create_mask(getfield(struct,field_name));
       case {'Norm'}
         im=tom_norm(im,getfield(getfield(struct,field_name),'Value'));
       case {'Cut'}
         val_tmp=getfield(struct,field_name),'Value') 
         im=tom_cut_out(im,val_tmp(1:(size(val_tmp,2)./2)),val_tmp((size(val_tmp,2)./2):(size(val_tmp,2)) );
       case {'Bin'}
           im=tom_bin(im,getfield(getfield(struct,field_name),'Value'));
       otherwise
           error([getfield(struct,field_name) ' corrupted']);
   end;

  
 end;
 
  if (struct.stack==1)
        out(:,:,ii)=im;  
   end
end;


