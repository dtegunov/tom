function tom_av2_eman_align2box(f_align,box_folder,box_size,funny_number)
%TOM_AV2_EMAN_ALIGN2BOX transfers a tom align2d file 2 a eman box file 
%
%   tom_av2_eman_align2box(f_align,box_folder,box_size,funny_number)
%
%  TOM_AV2_EMAN_ALIGN2BOX transfers a tom align2d file 2 a eman 
%  box file  
%   
%
%PARAMETERS
%
%  INPUT
%   f_align        filename of the alignment file
%   box_folder     foldername where all the box files will be written
%   box_size       (align2d.radius) vector
%   funny_number   (-3) addinal number in box file
%
%                        
%
%  EXAMPLE
%     
%  tom_av2_eman_align2box('picklist_high.mat','box_files_low');
%   
%  tom_av2_eman_align2box('picklist_high.mat','box_files_low',[360 360],-3);
%   
%   
%
%REFERENCES
%
%SEE ALSO
%   tom_aling2d
%
%   created  by fb
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



if (nargin<4)
    funny_number=-3;
end;

%load and initialize !
load(f_align);


if (isempty(align2d(1,1).radius) && exist('box_size','var')==0)
    errordlg('no radius or box size defined!!!');
end;

if (exist('box_size','var'))
    align2d(1,1).radius=round(box_size(1)./2);
else
    box_size=[align2d(1,1).radius.*2 align2d(1,1).radius.*2];
end;

%get unique list of filenames
for i=1:size(align2d,2)
    all_names{i}=align2d(1,i).filename;
end;
name_list=unique(all_names);


zz=1;
for i=1:length(name_list)
    idx=find(ismember(all_names,name_list{i}));
    for ii=1:length(idx)
        coords_x(ii)=align2d(1,idx(ii)).position.x-round(align2d(1,1).radius);
        coords_y(ii)=align2d(1,idx(ii)).position.y-round(align2d(1,1).radius);
    end;
    [a b c]=fileparts(name_list{i});
    write_box([box_folder '/' b '.box'],coords_x,coords_y,box_size,funny_number);
    disp([name_list{i} ' ' num2str(zz) ' of  ' num2str(length(name_list)) ]);
    clear('coords_x');
    clear('coords_y');
    zz=zz+1;
end;




function write_box(filename,coords_x,coords_y,box_size,funny_number)

fid=fopen(filename,'w');

for i=1:length(coords_x)
    fprintf(fid,'%d\t%d\t%d\t%d\t%d\n',coords_x(i),coords_y(i),box_size(1),box_size(2),funny_number);
end;

fclose(fid);