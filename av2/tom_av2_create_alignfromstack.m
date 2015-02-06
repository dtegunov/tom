function align2d = tom_av2_create_alignfromstack(stackfilename,outputname)
%TOM_AV2_CREATE_ALIGNFROMSTACK creates a particle alignment file from a stack
%
%   align = tom_av2_create_alignfromstack(stackfilename)
%
%   Note: the output alignment file is only suitable for
%    tom_av2_stackbrowser, not usable to generate a new stack with
%    tom_av2_particlepickergui.
%
%PARAMETERS
%
%  INPUT
%   stackfilename       full path and filename of the input stack or size vect
%   outputname          name for hd output   
%  
%  OUTPUT
%   align               alignment structure with default values
%
%EXAMPLE
%   align2d = tom_av2_create_alignfromstack('class_3.em');
%   save class_3.mat align2d
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by AK 09/02/06
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

if (nargin<2)
    outputname='';
end;


if (ischar(stackfilename))
    stackheader = tom_reademheader(stackfilename);
    if (isnumeric(outputname))
        [a b c]=fileparts(stackfilename);
        if (isempty(a))
            a='.';
        end;
        outputname=[a '/' b '.mat'];
    end;
else
    stackheader.Header.Size=stackfilename;
end;


for pointnumber=1:stackheader.Header.Size(3)
    align2d(1,pointnumber).dataset = '';
    align2d(1,pointnumber).filename = '';
    align2d(1,pointnumber).position.x = 0;
    align2d(1,pointnumber).position.y = 0;
    align2d(1,pointnumber).class = 'default';
    align2d(1,pointnumber).radius = stackheader.Header.Size(1)./2;
    align2d(1,pointnumber).color = [0 1 0];
    align2d(1,pointnumber).shift.x = 0;
    align2d(1,pointnumber).shift.y = 0;
    align2d(1,pointnumber).angle = 0;
    align2d(1,pointnumber).isaligned = 0;
    align2d(1,pointnumber).ccc = 0;
    align2d(1,pointnumber).quality = 0;
    align2d(1,pointnumber).normed = 'none';
    align2d(1,pointnumber).ref_class = 0;
end


if (isempty(outputname)==0)
    save(outputname,'align2d');
end;
