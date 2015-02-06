function out=tom_codonics_overview(filename,range1,range2,number,use)
%TOM_CODONICS_OVERVIEW prints a matrix of images to find the right parameters
%
%   out=tom_codonics_overview(filename,range1,range2,number,use)
%
%PARAMETERS
%
%  INPUT
%   filename            ...
%   range1              ...
%   range2              ...
%   number              ...
%   use                 ...
%  
%  OUTPUT
%   out                 ...
%
%EXAMPLE
%   tom_amira_createisosurface(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   tom_codonics_gui, tom_codonics
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


info=imfinfo(filename);

im_width=(info.Width+10).*number(1);
im_height=(info.Height+10).*number(2);

diff_width=im_width-2400;
diff_height=im_height-2680;

if (diff_width > diff_height)
    diff_max=diff_width;
    length=im_width;    
    length_max=2400;
else
    diff_max=diff_height;
    length=im_height;  
    length_max=2680;
end;

if (diff_max >= 0)
    scaleing=(length_max-2)./length;
    im_width=2400-1;
    im_height=2680-1;
else    
    scaleing=1;
end;


zi=1;
if (use(1)==1)
    if (zi==1)
        var_par1= ' CONTRAST ';
    else
        var_par2= ' CONTRAST ';
    end;
    zi=zi+1;
end;

if (use(2)==1)
    if (zi==1)
        var_par1=' GAMMA ';
    else
        var_par2=' GAMMA ';
    end;
    zi=zi+1;
end;

if (use(3)==1)
    if (zi==1)
        var_par1=' TCR ';
    else
        var_par2=' TCR ';
    end;
    
    zi=zi+1;
end;

if (use(4)==1)
    if (zi==1)
        var_par1=' MCM ';
    else
        var_par2=' MCM ';
    end;
    zi=zi+1;
end;

canvas=['CANVAS ' num2str(im_width) ' ' num2str(im_height) ];

fid=fopen('tmpppp1.txt','w');
fprintf(fid,'%c',canvas);
fclose(fid);


offset_x=0;
offset_y=0;
val1=range1(1);
val1_inkre=(range1(2)-range1(1))./(number(1)-1);
val2=range2(1);
val2_inkre=(range2(2)-range2(1))./((number(2)-1)) ;
zz=1;

for k=1:number(2)
    for i=1:number(1)
        place=['PLACE ' num2str(offset_x) ' ' num2str(offset_y)  var_par1 num2str(val1) var_par2 num2str(val2)...
                ' SCALE ' num2str(scaleing) ' BILINEAR '];
        out(k,i,1)=val1;
        out(k,i,2)=val2;
        fname=['tmppp_place_' num2str(zz) '.txt'];
        fid=fopen(fname,'w');
        fprintf(fid,'%c',place);
        fclose(fid);
        offset_x=offset_x+floor((info.Width+10).*scaleing);    
        val1=val1+val1_inkre;    
        zz=zz+1;    
    end;
    val2=val2+val2_inkre;
    val1=range1(1);
    offset_x=0;
    offset_y=offset_y+floor((info.Height+10).*scaleing);
end;


fid=fopen('tmpppp2.txt','w');
fprintf(fid,'PRINT');
fclose(fid);

fid=fopen('ftp_tmpppp128.txt','w');
% username
fprintf(fid,'nickell\n');
% 10 for free format print
fprintf(fid,'10\n');
fprintf(fid,'BINARY\n');
fprintf(fid,'put tmpppp1.txt\n');
for i=1:(number(1).*number(2))
    tmp_name=['put tmppp_place_' num2str(i) '.txt\n'];
    fprintf(fid,tmp_name);
    fprintf(fid,'put %s\n',filename');
end;

fprintf(fid,'put tmpppp2.txt\n');
fprintf(fid,'bye\n');
fclose(fid);


[s,w] = system('ftp -s:ftp_tmpppp128.txt codonics');

system('del tmpppp1.txt');
system('del tmpppp2.txt');

for i=1:(number(1).*number(2))
    tmp_name=['del tmppp_place_' num2str(i) '.txt'];
    system(tmp_name);
end;
system('del ftp_tmpppp128.txt');