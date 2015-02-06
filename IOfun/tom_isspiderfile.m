function out = tom_isspiderfile(filename)
%TOM_ISSPIDERFILE creates ...
%
%   out = tom_isspiderfile(filename)
%
%PARAMETERS
%
%  INPUT
%   filename            ...
%  
%  OUTPUT
%   out                 ...
%
%EXAMPLE
%   tom_isspiderfile(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   TOM_SPIDERWRITE, TOM_SPIDERHEADER, TOM_SPIDERREAD
%
%   created by AK 04/25/06
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


fid = fopen(filename,'rb','ieee-be');
if fid==-1
    out=-1;
    return;
end

test = fread(fid,5,'float');

if (isempty(test)==1)
    out=-1;
    fclose(fid);
    return;
end;
    
try
    fid = fopen(filename,'rb','ieee-le');
    testval = fread(fid,1,'float');
    fseek(fid,0,-1);
    if testval < 1 || testval > 10000
        fclose(fid);
        fid = fopen(filename,'rb','ieee-be');
        testval = fread(fid,1,'float');
        fseek(fid,0,-1);
        if testval < 1 || testval > 10000
            out=0;
        else
            out=1;
        end;        
    else
        out=1;
    end;
    fclose(fid);
    
catch
    error(['Could not open ' filename]);
    fclose(fid);
end




