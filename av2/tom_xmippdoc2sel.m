function tom_xmippdoc2sel(doc,sel)

% TOM_XMIPPDOC2SEL extracts the particles from the Xmipp-Doc file
% and creates a sel file.
%
% Input:
%           doc:                        Xmipp doc file.
% Output:
%           sel:                        Xmipp sel file.
% 
% Example:
% tom_xmippdoc2sel('all_parts.doc','all_parts.sel')
%
%SEE ALSO
%   tom_xmippdocread
%
%   created by SN 05/06/10
%
%   Nickell et al., 'TOM software toolbox: acquisition and analysis for electron tomography',
%   Journal of Structural Biology, 149 (2005), 227-234.
%
%   Copyright (c) 2004-2010
%   TOM toolbox for Electron Tomography
%   Max-Planck-Institute of Biochemistry
%   Dept. Molecular Structural Biology
%   82152 Martinsried, Germany
%   http://www.biochem.mpg.de/tom


fid_doc=fopen(doc);
fid_sel=fopen(sel,'w');

Header=fgetl(fid_doc);

while 1
    part=fgetl(fid_doc);
    if ~ischar(part), break, end
    align=fgetl(fid_doc);    
    fprintf(fid_sel,'%s %i\n',part(4:end),1);
end;
fclose(fid_doc);
fclose(fid_sel);
