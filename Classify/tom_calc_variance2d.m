function [var_st mea_st]=tom_calc_variance2d(stack,mask,norm,flag)
%tom_calc_variance calculates the variance of a feature vector
%
%    var=tom_calc_variance(input_data,classes,class)
%
%   plots positions as a 2d grid
%   useful for ploting a som
%
%PARAMETERS
%
%  INPUT
%   input_data    2d part stack
%   mask          mask or 'default' for default mask (rad-1), no_mask        
%   norm_flg      mean0+1std,phase... or 'default' (mean0+1std), no_norm
%                 for no norm  
%   flag          eco (save memory) or fast (calc paralell)
%  
%  OUTPUT
%   var           variance image
%
%EXAMPLE
%    var=tom_calc_variance(stack);
%
%REFERENCES
%
%SEE ALSO
%   
%
%   created by FB 11/17/09
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


sz=size(stack);

if (length(sz)==2)
    sz(3)=1;
end;

if (nargin <2);
    mask='default';
end;


if (nargin <3 || strcmp(norm,'default'));
    norm='mean0+1std';
end;


if (nargin <4);
    flag='eco';
end;

if (strcmp(mask,'default'))
    mask=tom_spheremask(ones(sz(1),sz(2)),round(sz(1)./2)-1 );
end;

if (strcmp(mask,'no_mask'))
    mask=ones(sz,sz);
end;
    

if (strcmp(norm,'no_norm')==0)
    parfor i=1:sz(3)
        stack(:,:,i)=tom_norm(stack(:,:,i),norm,mask).*mask;
    end;
end;

mea_st=sum(stack,3)./sz(3);

if (strcmp(flag,'eco'))
    var_st=zeros(sz(1),sz(2));
    for i=1:size(stack,3)
        var_st=var_st+(stack(:,:,i)-mea_st).^2;
    end;
else
    std_stack=zeros(sz);
    parfor i=1:sz(3)
        std_stack(:,:,i)=(stack(:,:,i)-mea_st).^2;
    end;
    var_st=sum(std_stack,3);
end;

if ((sz(3)-1)>0)
    var_st=var_st./(sz(3)-1);
else
   var_st=var_st; 
end;


if (strcmp(norm,'no_norm')==0)
     var_st=tom_norm(var_st,norm,mask).*mask;
end;







