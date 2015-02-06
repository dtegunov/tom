function out=tom_remove_nan(in,mean)
%TOM_REMOVE_NAN creates ...
%
%   out=tom_remove_nan(in,mean)
%
%PARAMETERS
%
%  INPUT
%   in                  ...
%   mean                ...
%  
%  OUTPUT
%   out                 ...
%
%EXAMPLE
%   ... = tom_remove_nan(...);
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


out=zeros(size(in));
st=zeros(size(in,3).*size(in,2).*size(in,1),1);

if nargin==1
    zz=1;
    for z=1:size(in,3)
        for y=1:size(in,2)
            for x=1:size(in,1)

                if (isnan(in(x,y,z))==0)
                    st(zz)=in(x,y,z);
                    zz=zz+1;
                else
                   
                end;
            end;
        end;
    
    end;
    m=tom_dev(st,'noinfo');
else
    m=mean;
end;



for z=1:size(in,3)
    for y=1:size(in,2)
        for x=1:size(in,1)
            
            if (isnan(in(x,y,z))==0)
                out(x,y,z)=in(x,y,z);
            else
                out(x,y,z)=m;
            end;
        
        end;
    end;

end;