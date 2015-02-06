function stack=tom_file2stack(path,num,ext,ibin)
%TOM_FILE2STACK creates ...
%
%   stack=tom_file2stack(path,num,ext)
%
%PARAMETERS
%
%  INPUT
%   path                name of files that are converted (without index and
%                       extension
%   num                 max index of files to be stored
%   ext                 extension of input files
%   ibin                binning (optional)
%  
%  OUTPUT
%   stack               ...
%
%EXAMPLE
%   ... = tom_file2stack(...);
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
if nargin < 4
    ibin=0;
end;

for i=1:num
    try
        im=tom_emreadc([path num2str(i) ext]);
        if ibin>0
            stack(:,:,i)=tom_bin(im.Value,ibin);
        else
            stack(:,:,i)=im.Value;
        end;
    catch
        disp([path num2str(i) ext ' not found']);
    end;
    
    
end;