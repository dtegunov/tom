function out=tom_iseman2h5file(filename)
% tom_iseman2h5file tests for a eman2,sparx hdf5 file
%
% [out] = tom_isdm3file(filename)
%
% PARAMETERS
% 
%   INPUT
%    filename                filename of the file
%
%   OUTPUT
%    out                     1 is e2,sparx 
%                            0 is not e2,sparx
%
%SEE ALSO
%     TOM_ISEMFILE,TOM_EMAN2_READ
%
%   created by fb 06/09/12
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

out=1;
loc_num_of_entry='/MDF/images/';
att_name_num_of_entry='imageid_max';

try
    num_of_entries = h5readatt(filename,loc_num_of_entry,att_name_num_of_entry);
catch Me
    out=0;
end;



