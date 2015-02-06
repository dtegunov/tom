function [out] = tom_isxmipp_fscfile(em_name)
%  tom_isxmipp_fscfile tests for a file in Xmipp FSC-Format
%
%  [out] = tom_isxmipp_fscfile(em_name)
%
%    Checks if the file is a Xmipp FSC-Format,
%	 a text format with a Fourier Shell Correlation out of Xmipp.
%
%    Structure of Xmipp_FSC-Data Files:
%   # Resol. [1/Ang]      FRC      FRC_random_noise     Resol. [Ang]
%   0.00176753           0.9997           0.5547           565.76
%   0.00353507           0.999702         0.328798         282.88
%   0.0053026            0.99936          0.264906         188.587 ...
%
%  INPUT
%   em_name             Filename
%
%  OUTPUT
%   out                 is 1 for true, yes the file is in EM-format.
%                       is 0 for false, no the file is not in EM-format.
%                       is -1 for file doesn't exist at all.
%
%SEE ALSO
%     TOM_EMWRITE, TOM_EMHEADER, TOM_EMREADEMHEADER
%
%   created by SN 07/23/09
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


% open the stream with the correct format !
fid = fopen(em_name,'r');
if fid==-1
    out=-1;
    return;
end;
C=textscan(fid,'%c',2);
fclose(fid);
if (isempty(C))
    out=-1;
    return;
end;
c=C{1};
if findstr(c','#R')
    out=1;
else
    out=0;
end;

    



