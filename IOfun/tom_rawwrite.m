function tom_rawwrite(em_name,Data,endian)
% TOM_RAWWRITE writes data in an raw-file format
%
%   tom_rawwrite(em_name,Data,endian)
%
%PARAMETERS
%
%  INPUT
%   name                ['PATHNAME' 'FILENAME']
%   Data                Values of Image Data
%   endian              Endianess, PC ... little endian 'l',
%                       big endian 'b'
%  OUTPUT
%
%EXAMPLE
%
%   tom_rawwrite('HPIEMV.raw',out);
%
%
%REFERENCES
%
%SEE ALSO
%   TOM_RAWREAD,TOM_EMWRITE, TOM_EMHEADER, TOM_READEMHEADER
%
%   created by SN 12/01/09
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


if nargin <1 error(['Data not specified']);  end;

if nargin <3 

    [~,~,endian] = computer;

end;

fid = fopen(em_name,'w',endian); 

count=fwrite(fid,Data,class(Data));

fclose(fid);

if count ~= numel(Data)
    error('Cannot write file in tom_rawwrite.');
end;

