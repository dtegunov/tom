function outlist = tom_av3_align_sum(Align,flag)
%TOM_AV3_ALIGN_SUM sums the shifts and angles of a align list, needed for averaging
%
%   outlist = tom_av3_align_sum(Align)
%
%PARAMETERS
%
%  INPUT
%   Align               ...
%  
%  OUTPUT
%   outlist             ...
%
%EXAMPLE
%   ... = tom_av3_align_sum(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   TOM_AV3_PASTE4OSCARGUI, TOM_AV3_OSCARGUI, TOM_AV_EXTRACT_ANGLESSHIFTS,
%   TOM_AV_EXTRACTANGLESSHIFTSGUI, TOM_AV3_AVERAGE
%
%   created by AK 11/10/05
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
    flag='trans_rot';
end;

outlist = Align(end,:);

%loop over particles
for i = 1:size(Align,2)
    %loop over runs
    for ii=1:size(Align,1);
        shifts_in(ii,1)=Align(ii,i).Shift.X;
        shifts_in(ii,2)=Align(ii,i).Shift.Y;
        shifts_in(ii,3)=Align(ii,i).Shift.Z;
        euler_in(ii,1)=Align(ii,i).Angle.Phi;
        euler_in(ii,2)=Align(ii,i).Angle.Psi;
        euler_in(ii,3)=Align(ii,i).Angle.Theta;
    end;
    
    [euler_out,shifts_out,rotmatrix]=tom_sum_rotation(euler_in,shifts_in,flag);
    outlist(1,i).Shift.X = shifts_out(1);
    outlist(1,i).Shift.Y = shifts_out(2);
    outlist(1,i).Shift.Z = shifts_out(3);
    %rotation matrix
    outlist(1,i).Angle.Rotmatrix = rotmatrix;
    outlist(1,i).Angle.Phi = euler_out(1);    
    outlist(1,i).Angle.Psi = euler_out(2);
    outlist(1,i).Angle.Theta = euler_out(3);
end