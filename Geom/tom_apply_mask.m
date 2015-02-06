function im=tom_apply_mask(im,mask,text,flag)
%TOM_APPLY_MASK performs exception handling for tom_av2_rec gui
%
%   im=tom_apply_mask(im,mask,text,flag)
%
%PARAMETERS
%
%  INPUT
%   im                  handles struct of gui
%   mask                struct including gui values
%   text                string with items to be checked
%   flag                error or warning text
%  
%  OUTPUT
%   im                  image.*mask
%
%DESCRIPTION
%   happy young sporty function searching a nice gui for  
%   potential errorhandling
%
%EXAMPLE
%   im=tom_emread('pyrodictium_1.em');
%   im=tom_apply_mask(im,tom_spheremask(ones(30,30),5),'interactive-save','mas
%k for projection differs in size');
%
%REFERENCES
%
%SEE ALSO
%   get_gui_values
%
%   created by FB 20/10/06
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


if (nargin < 2)
    error('at least 2 parameters');
end;
    
if (nargin < 3)
    text='mask and image differ in size';
end;

if (nargin < 4)
    flag='save-interactive';
end;



if (strcmp(flag,'save-interactive'))
    dims_im=max(size(size(im)));
    dims_mask=max(size(size(im)));
    
    if (dims_im~=dims_mask)
        im=0;
        errordlg(text);
    end;
    
    if (sum(size(mask)==size(im))~=dims_im)
        im=0;
        errordlg(text);
    end;
    im=im.*mask;

end;


if (strcmp(flag,'wild-west'))
    
end;

