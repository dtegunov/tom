function [matlab_surface matlab_normals]=tom_stl2surface(stl_filename)

% TOM_STL2SURFACE reads a .stl file into a matlab surface structure.
%
%   [matlab_surface matlab_normals]=tom_stl2surface(stl_filename)
%
% PARAMETERS
%   INPUT
%   stl_filename    name of the .stl file 
%
%   OUTPUT
%   matlab_surface  surface
%   matlab_normals  surface normals
%
% EXAMPLE
%   % You can create the simplex.stl file with the tom_surface2stl example
%   [surface normals]=tom_stl2surface('simplex.stl');
%   p=patch('Faces',surface.faces,'Vertices',surface.vertices,'FaceVertexCData'...
%   ,ones(size(surface.faces,1),3).*.8,'FaceColor','flat');
%   view(3);grid on; axis image;
%   set(p,'FaceColor', [230./255 160./255 45./255],'EdgeColor', ...
%   [230./255 160./255 45./255])
%   material shiny; light;
%   lighting phong; set(p,'EdgeLighting','phong')
%
%REFERENCES
%
% SEE ALSO
%
%   last change
%   07/20/08 SN
%
%   Nickell et al., 'TOM software toolbox: acquisition and analysis for electron tomography',
%   Journal of Structural Biology, 149 (2005), 227-234.
%
%   Copyright (c) 2004-2008
%   TOM toolbox for Electron Tomography
%   Max-Planck-Institute of Biochemistry
%   Dept. Molecular Structural Biology
%   82152 Martinsried, Germany
%   http://www.biochem.mpg.de/tom


fid = fopen(stl_filename);
if fid==-1
    error('cannot open .stl file.'); 
end;
s=fread(fid,5,'char');
fseek(fid,0,-1);
if isequal(char(s'),'solid') % then its ASCII format

    % here is the ASCII format
    [matlab_vertex,matlab_normals,fid]=read_stl_file_ascii(fid);

else % here is binary format
    [nr_faces,fid]=read_header_stl_file_binary(fid);
    [matlab_vertex,matlab_normals,fid]=read_stl_file_binary(nr_faces, fid);
    matlab_vertex=matlab_vertex';
end;
fclose(fid);
[matlab_surface.vertices b c]=unique(matlab_vertex','rows');;
matlab_surface.faces=reshape(c,[3 size(c,1)./3])';


function found=find_normal(fid)
in=''; found=0;
while ~isequal(in,'normal') && ~feof(fid)
    in=fscanf(fid,'%s',1);
    found=1;
end;
function find_vertex(fid)
in='';
while ~isequal(in,'vertex') && ~feof(fid)
    in=fscanf(fid,'%s',1);
end;


function [nr_faces,fid]=read_header_stl_file_binary(fid)
tmp=fread(fid,80,'uchar');
nr_faces=fread(fid,1,'unsigned long');

function [vertices,normals,fid]=read_stl_file_binary(nr_faces, fid);
normals=zeros(nr_faces,3);
vertices=zeros(nr_faces.*3,3);
iii=1;
for i=1:nr_faces
    normals(i,:)=fread(fid,3,'float');
    for ii=1:3;
        vertices(iii,:)=fread(fid,3,'float');
        iii=iii+1;
    end;
    fread(fid,1,'short');
end;

function [vertices,normals,fid]=read_stl_file_ascii(fid);

normals=[]; vertices=[];

while ~feof(fid) && find_normal(fid)
    normal=fscanf(fid,'%f %f %f');
    if isempty(normal); break; end;
    normals=[normals, normal];
    for i=1:3
        find_vertex(fid);
        vertex=fscanf(fid,'%f %f %f');
        vertices=[vertex, vertices];
    end;
end;
