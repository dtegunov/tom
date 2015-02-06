function img_out=tom_reshape_stack_memory(stack,binning,mask,normstr)
%TOM_RESHAPE_STACK_MEMORY reshapes stack
%
%   img_out=tom_reshape_stack_memory(stack,binning,mask,normstr);
%PARAMETERS
%
%  INPUT
%   stack               stack of particles
%   binning             binning
%   mask                mask
%   normstr             normstr 'mean0+1std' '3std' ...
%
%  
%  OUTPUT
%   img_out             reshaped stack
%
%EXAMPLE
%   img_out=tom_reshape_stack_memory(stack,1);
%   
%  
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by ... fb(eckster) anno 2008
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

if nargin < 2
    binning=0;
end;

if nargin < 3 
    mask=ones(size(stack,1),size(stack,2));
end;

if nargin < 4
    normstr='no_norm';
end;


sz=size(stack);

if (binning >0)
    sz(1)=floor(sz(1)./2^binning);
    sz(2)=floor(sz(2)./2^binning);
end;

img_out=zeros(sz(1).*sz(2),sz(3));

for i=1:sz(3)
    tmp=tom_bin(stack(:,:,i),binning).*tom_bin(mask,binning);
    if strcmp(normstr,'no_norm')==0
        tmp=tom_norm((tmp+2).*2,normstr,tom_bin(mask,binning));
    end;
    img_out(:,i)=reshape(tmp,sz(1).*sz(2),'');
end;