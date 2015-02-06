function tom_mrc2em(pathname, newname)
%TOM_MRC2EM converts FEI style MRC-format files to EM format
%
%   tom_mrc2em(pathname, newname)
%
%PARAMETERS
%
%  INPUT
%   pathname            ...
%   newname             ...
%  
%  OUTPUT
%
%DESCRIPTION
%   tom_mrc2em opens a dialog box where a folder can be selected. This folder 
%   should contains MRC file only. All the MRC files will be converted to EM format.
%   This should be used to convert a full MRC tilt series.
%
%   tom_mrc2em(name_mrc_file, name_em_file) converts name_mrc_file to name_em_file.
%   name_mrc_file contains the path and name of the mrc file. 
%   name_em_file contains the path and name of the destination file. 
%   Can convert 2D or 3D image.
%
%EXAMPLE
%   tom_mrc2em  a dialog box appears. Select the folder containing the mrc files
%
%   tom_mrc2em('c:\test\mrs_001.mrc', 'c:\test\mrs_001.em') converts MRC file mrs_001.mrc to mrs_001.em (EM Format)
%
%REFERENCES
%
%SEE ALSO
%   TOM_MRCREAD, TOM_MRCSTACK2EMSERIES, TOM_EMWRITE
%
%   created by SN 08/07/02
%   updated by WDN 05/23/05
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


if nargin <1 
    path_org=pwd;
    [pathname] = uigetdir(path_org,'Select the source-directory of your tiltseries');
    if isequal(pathname,0) 
        disp('Cancelled.');
        return;
    end
    newextension='.em';
    pathname=[pathname '\'];
    d=dir(pathname);
    laufx=0;
    pathname
    for lauf=1:(size(d,1)-2)
        mrcname=d(lauf+2).name;
        sfn=size(mrcname,2);
        if strcmp(mrcname(sfn-2:sfn), 'mrc')
            newemname=[mrcname(1:sfn-3) 'em'];
        else
            newemname=[mrcname '.em'];
        end
        i=tom_mrcread([pathname d(lauf+2).name],'le');
        n=i.Value;
        n=tom_emheader(n);
        n.Header.Tiltangle=i.Header.Tiltangle;
        n.Header.Voltage=300000;
        n.Header.Cs=2.0;
        n.Header.Aperture=0;
        n.Header.Magnification=i.Header.Magnification;
        %n.Header.Postmagnification=i.Header.Postmagnification;
        n.Header.Tiltaxis=i.Header.Tiltaxis;
        n.Header.Exposuretime=i.Header.Exposuretime;
        n.Header.Objectpixelsize=i.Header.Objectpixelsize*1e9;
        n.Header.Microscope='Polara         ';
        n.Header.Defocus=i.Header.Defocus;
        n.Header.Pixelsize=30000;
        n.Header.CCDArea=n.Header.Pixelsize*i.Header.Size(1);
        tom_emwrite([pathname newemname],n);
        disp([d(lauf+2).name ' --> ' newemname ' ']);
    end;
    disp('Conversion done');   
end

if nargin==1
    error('Not enough input parameter');
    n=0;i=0;
end

if nargin==2
    i=tom_mrcread(pathname,'le');
    n=i.Value;
    n=tom_emheader(n);
    n.Header.Tiltangle=i.Header.Tiltangle;
    n.Header.Voltage=300000;
    n.Header.Cs=2.0;
    n.Header.Aperture=0;
    n.Header.Magnification=i.Header.Magnification;
    %n.Header.Postmagnification=i.Header.Postmagnification;
    n.Header.Tiltaxis=i.Header.Tiltaxis;
    n.Header.Exposuretime=i.Header.Exposuretime;
    n.Header.Objectpixelsize=i.Header.Objectpixelsize*1e9;
    n.Header.Microscope='Polara         ';
    n.Header.Defocus=i.Header.Defocus;
    n.Header.Pixelsize=30000;
    n.Header.CCDArea=n.Header.Pixelsize*i.Header.Size(1);
    if (isempty(findstr(newname,'.em'))) & (isempty(findstr(newname,'.vol')))
        if size(size(n.Value),2)<=2;
            newname=[newname '.em'];
        else size(size(n.Value),2)==3;
            newname=[newname '.vol'];
        end       
    end
    tom_emwrite([newname],n);
    disp([pathname ' >> ' newname ' ']);
end
clear i;
clear n;

