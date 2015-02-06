function tom_codonics(filename,gamma,contrast,tcr,mcm,size,default)
%TOM_CODONICS prints an image
%
%   tom_codonics(filename,gamma,contrast,tcr,mcm,size,default)
%
%   tom_codonics prints an image using ftp control.
%
%PARAMETERS
%
%  INPUT
%   filename            name of the tif image which should be printed
%   gamma               gamma value (darker or lighter image) (see also codonics manual)
%   contrast            contrast value (see also codonics manual)
%   tcr                 tcr value true color randering (see also codonics manual)
%   mcm                 mcm value medical color matching (see also codonics manual)
%   size                vector determining the size of the image in cm [widht hight]
%   default             vector including the information whether to
%                        use the given values for (gamma, contrast,
%                        tcr mcm,size) or to use the default value. 0
%                        for using default
%  
%  OUTPUT
%
%EXAMPLE
%   tom_codonics('c:\test',3,1,5,2,[5 5],[1 0 1 1 1])
%
%REFERENCES
%
%SEE ALSO
%   tom_codonics_gui, tom_codonics_overview
%
%   created by FB 04/03/04
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

if (nargin==0)
    tom_codonics_gui;
    return;
end;


info=imfinfo(filename);

width_dest=round(118.*size(1));
hight_dest=round(118.*size(2));
width_is=info.Width;
hight_is=info.Height;

im=imread(filename);

if (default(1)==0)
    im=imresize(im,[hight_dest width_dest],'bilinear');
end;
imwrite(im,'im_tmp.tif','tiff');

if (default(1)==0)
    canvas=['CANVAS ' num2str(width_dest+1) ' ' num2str(hight_dest+1)];
    
else
    canvas=['CANVAS ' num2str(width_is+1) ' ' num2str(hight_is+1)];
    
end;

fid=fopen('tmpppp1.txt','w');
fprintf(fid,'%c',canvas);
fclose(fid);
place=['PLACE 0 0'];

if (default(2)==0)
    place=[place ' CONTRAST '  num2str(contrast)];
end;

if (default(3)==0)
    place=[place ' GAMMA ' num2str(gamma) ]; 
end;

if (default(4)==0)
    place=[place ' TCR ' num2str(tcr) ]; 
end;

if (default(5)==0)
    place=[place ' MCM ' num2str(mcm) ]; 
end;


fid=fopen('tmpppp2.txt','w');
fprintf(fid,'%c',place);
fclose(fid);
fid=fopen('tmpppp3.txt','w');
fprintf(fid,'PRINT');
fclose(fid);

fid=fopen('ftp_tmpppp128.txt','w');
% username
fprintf(fid,'nickell\n');
% 10 for free format print
fprintf(fid,'10\n');
fprintf(fid,'BINARY\n');
fprintf(fid,'put tmpppp1.txt\n');
fprintf(fid,'put tmpppp2.txt\n');
fprintf(fid,'put %s\n','im_tmp.tif');
fprintf(fid,'put tmpppp3.txt\n');
fprintf(fid,'bye\n');
fclose(fid);


[s,w] = system('ftp -s:ftp_tmpppp128.txt codonics');

delete('tmppp*.txt');
delete('im_tmp.tif');
delete('ftp_tmpppp128.txt');