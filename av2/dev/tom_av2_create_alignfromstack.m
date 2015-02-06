function align = tom_av2_create_alignfromstack(stackfilename)
% creates a particle alignment file from a stack
%
% Note: the output alignment file is only suitable for
% tom_av2_stackbrowser, not usable to generate a new stack with
% tom_av2_particlepickergui.
%
% SYNTAX
% align = tom_av2_create_alignfromstack(stackfilename)
%
% INPUT
% stackfilename:    full path and filename of the input stack
%
% OUTPUT
% align:            alignment structure with default values
%
%Copyright (c) 2006
%TOM toolbox for Electron Tomography
%Max-Planck-Institute for Biochemistry
%Dept. Molecular Structural Biology
%82152 Martinsried, Germany
%http://www.biochem.mpg.de/tom
%
%Created: 09/02/06 AK

stackheader = tom_reademheader(stackfilename);

for pointnumber=1:stackheader.Header.Size(3)
    align(1,pointnumber).dataset = '';
    align(1,pointnumber).filename = '';
    align(1,pointnumber).position.x = 0;
    align(1,pointnumber).position.y = 0;
    align(1,pointnumber).class = 'default';
    align(1,pointnumber).radius = stackheader.Header.Size(1)./2;
    align(1,pointnumber).color = [0 1 0];
    align(1,pointnumber).shift.x = 0;
    align(1,pointnumber).shift.y = 0;
    align(1,pointnumber).angle = 0;
    align(1,pointnumber).isaligned = 0;
    align(1,pointnumber).ccc = 0;
    align(1,pointnumber).quality = 0;
    align(1,pointnumber).normed = 'none';
end