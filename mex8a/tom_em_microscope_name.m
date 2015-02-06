function res = tom_em_microscope_name(mic)
%TOM_EM_MICROSCOPENAME Converts the number in the EM-Header representing
%  the microscope to a string and back
%
%res = tom_em_microscope_name(mic)
%
%EXAMPLE
%  name = tom_em_microscope_name(7);
%  number = tom_em_microscope_name('Titan          ');
%
%SEE ALSO
%   TOM_EMREADC3
%
%   created by Thomas Haller Nov. 13 2007
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



name_list = { 'EM420          ', ...
              'CM12           ', ...
              'CM200          ', ...
              'CM120/Biofilter', ...
              'CM300          ', ...
              'Polara         ', ...
              'Titan          ', ...
              'Tecnai F20     ', ...
              'extern         ' };

if (ischar(mic))
    res = find(strcmp(name_list, mic));
    if (isempty(res) || res == length(name_list))
        res = 0;
    end;
elseif (isnumeric(mic) && numel(mic)==1)
    if (mic ~= round(mic) || mic<=0 || mic>length(name_list))
        res = name_list{end};
    else
        res = name_list{mic};
    end;
else
    error('Give name or number.');
end;


