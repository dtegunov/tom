function [in mean_out]=tom_rm_mean(in,mean_in)
%TOM_RM_MEAN creates ...
%
%   [in mean_out]=tom_rm_mean(in,mean_in)
%
%PARAMETERS
%
%  INPUT
%   in                  ...
%   mean_in             ...
%  
%  OUTPUT
%   in          		...
%   mean_out      		...
%
%EXAMPLE
%   ... = tom_rm_mean(...);
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


if (nargin==1)

    parfor i=1:size(in,2)
        mean_out(i)=mean(in(:,i));
        in(:,i)=in(:,i)-mean_out(i);
    end;

else
    %inverse case
     for i=1:size(in,2)
        mean_out(i)=mean(in(:,i));
        in(:,i)=in(:,i)+mean_in(i);
    end;
    

end;