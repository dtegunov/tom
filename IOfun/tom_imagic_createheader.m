function out = tom_imagic_createheader(Data,type)
%TOM_IMAGIC_CREATEHEADER creates an imagic structure of a matrix
%
%   out = tom_imagic_createheader(Data,type)
%
%PARAMETERS
%
%  INPUT
%   Data                input matrix
%   type                must be '2dstack'
%  
%  OUTPUT
%   out                 imagic structure
%
%
%REFERENCES
%
%   http://www.imagescience.de/formats/formats.htm
%
%SEE ALSO
%   TOM_IMAGICREAD, TOM_IMAGICWRITE
%
%   created by AK 08/12/06
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

if isstruct(Data)
    error('Input Data is already a structure.');
end

out.Value = Data;
[y, m, d, h, mn, s] = datevec(datenum(now));

if isequal(lower(type),'2dstack')
    for i=1:size(Data,3)
        [a,b,c,d,e]=tom_dev(Data(:,:,i),'noinfo');
        out.Header.IMAGIC(i).IMN = i-1;
        out.Header.IMAGIC(i).IFOL = size(Data,3);
        out.Header.IMAGIC(i).IERROR = 0;
        out.Header.IMAGIC(i).NHFR = 1;
        out.Header.IMAGIC(i).NMONTH = m;
        out.Header.IMAGIC(i).NDAY = d;
        out.Header.IMAGIC(i).NYEAR = y;
        out.Header.IMAGIC(i).NHOUR = h;
        out.Header.IMAGIC(i).NMINUT = mn;
        out.Header.IMAGIC(i).NSEC = s;
        out.Header.IMAGIC(i).NPIX2 = size(Data,1).*size(Data,2);
        out.Header.IMAGIC(i).NPIXEL = size(Data,1).*size(Data,2);
        out.Header.IMAGIC(i).IXLP = size(Data,1);
        out.Header.IMAGIC(i).IYLP = size(Data,2);
        out.Header.IMAGIC(i).TYPE = 'REAL';
        out.Header.IMAGIC(i).IXOLD = 0;
        out.Header.IMAGIC(i).IYOLD = 0;
        out.Header.IMAGIC(i).AVDENS = a;
        out.Header.IMAGIC(i).SIGMA = d;
        out.Header.IMAGIC(i).VARIAN = e;
        out.Header.IMAGIC(i).OLDAVD = 0;
        out.Header.IMAGIC(i).DENSMAX = b;
        out.Header.IMAGIC(i).DENSMIN = c;
        out.Header.IMAGIC(i).COMPLEX = 0;
        out.Header.IMAGIC(i).CXLENGTH = 0;
        out.Header.IMAGIC(i).CYLENGTH = 0;
        out.Header.IMAGIC(i).CZLENGTH = 0;
        out.Header.IMAGIC(i).CALPHA = 0;
        out.Header.IMAGIC(i).CBETA = 0;
        out.Header.IMAGIC(i).NAME = '';
        out.Header.IMAGIC(i).CGAMMA = 0;
        out.Header.IMAGIC(i).MAPC = 0;
        out.Header.IMAGIC(i).MAPR = 0;
        out.Header.IMAGIC(i).MAPS = 0;
        out.Header.IMAGIC(i).ISPG = 0;
        out.Header.IMAGIC(i).NXSTART = 0;
        out.Header.IMAGIC(i).NYSTART = 0;
        out.Header.IMAGIC(i).NZSTART = 0;
        out.Header.IMAGIC(i).NXINTV = 0;
        out.Header.IMAGIC(i).NYINTV = 0;
        out.Header.IMAGIC(i).NZINTV = 0;
        out.Header.IMAGIC(i).IZLP = 0;
        out.Header.IMAGIC(i).I4LP = 0;
        out.Header.IMAGIC(i).I5LP = 0;
        out.Header.IMAGIC(i).I6LP = 0;
        out.Header.IMAGIC(i).ALPHA = 0;
        out.Header.IMAGIC(i).BETA = 0;
        out.Header.IMAGIC(i).GAMMA = 0;
        out.Header.IMAGIC(i).IMAVERS = 0;
        out.Header.IMAGIC(i).REALTYPE = 0;
        out.Header.IMAGIC(i).RONLY = 0;
        out.Header.IMAGIC(i).ANGLE = 0;
        out.Header.IMAGIC(i).RCP = 0;
        out.Header.IMAGIC(i).IXPEAK = 0;
        out.Header.IMAGIC(i).IYPEAK = 0;
        out.Header.IMAGIC(i).CCC = 0;
        out.Header.IMAGIC(i).ERRAR = 0;
        out.Header.IMAGIC(i).ERR3D = 0;
        out.Header.IMAGIC(i).REF = 0;
        out.Header.IMAGIC(i).CLASSNO = 0;
        out.Header.IMAGIC(i).LOCOLD = 0;
        out.Header.IMAGIC(i).OLDAVD = 0;
        out.Header.IMAGIC(i).OLDSIGMA = 0;
        out.Header.IMAGIC(i).XSHIFT = 0;
        out.Header.IMAGIC(i).YSHIFT = 0;
        out.Header.IMAGIC(i).NUMCLS = 0;
        out.Header.IMAGIC(i).OVQVAL = 0;
        out.Header.IMAGIC(i).EANGLE = 0;
        out.Header.IMAGIC(i).EXSHIFT = 0;
        out.Header.IMAGIC(i).EYSHIFT = 0;
        out.Header.IMAGIC(i).CMTOTVAR = 0;
        out.Header.IMAGIC(i).INFORMAT = 0;
        out.Header.IMAGIC(i).NUMEIGEN = 0;
        out.Header.IMAGIC(i).NIACTIVE = 0;
        out.Header.IMAGIC(i).RESOLX = 0;
        out.Header.IMAGIC(i).RESOLY = 0;
        out.Header.IMAGIC(i).RESOLZ = 0;
        out.Header.IMAGIC(i).ALPHA2 = 0;
        out.Header.IMAGIC(i).BETA2 = 0;
        out.Header.IMAGIC(i).GAMMA2 = 0;
        out.Header.IMAGIC(i).FABOSA1 = 0;
        out.Header.IMAGIC(i).FABOSA2 = 0;
        out.Header.IMAGIC(i).FABOSA3 = 0;
        out.Header.IMAGIC(i).NMETRIC = 0;
        out.Header.IMAGIC(i).ACTMSA = 0;
        out.Header.IMAGIC(i).COOSMSA = 0;
        out.Header.IMAGIC(i).EIGVAL = 0;
        out.Header.IMAGIC(i).HISTORY = '';
    end
end