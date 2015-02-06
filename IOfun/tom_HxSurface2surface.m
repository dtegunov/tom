function [matlab_surface]=tom_HxSurface2surface(hx_filename)

% TOM_HXSURFACE2SURFACE reads a HxSurface file (Amira, ASCII) into a matlab surface structure.
%
%   [matlab_surface]=tom_HxSurface2surface(hx_filename)
%
% PARAMETERS
%   INPUT
%   hx_filename    name of the HxSurface file 
%
%   OUTPUT
%   matlab_surface  surface
%
% EXAMPLE
%   % You can create the simplex.stl file with the tom_surface2HxSurface example
%   [surface]=tom_HxSurface2surface('simplex.surf');
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
%   tom_pov2dxf, tom_pov2dxf_gui, tom_surface2stl, tom_surface2dxf, patch,
%   reducepatch, tom_surface2HxSurface, tom_stl2surface
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


fid = fopen(hx_filename);
if fid==-1
    error('cannot open HxSurface file.'); 
end;
% read the vertices
find_vertex(fid);
in=fscanf(fid,'%s',1);
nr_vertices=str2num(in);
matlab_surface.vertices=zeros(nr_vertices,3);
for i=1:nr_vertices
    in=fscanf(fid,'%f',3);
    matlab_surface.vertices(i,:)=in;
end;
% read the faces
find_face(fid);
in=fscanf(fid,'%s',1);
nr_faces=str2num(in);
matlab_surface.faces=zeros(nr_faces,3);
for i=1:nr_faces
    in=fscanf(fid,'%f',3);
    matlab_surface.faces(i,:)=in;
end;
fclose(fid);

function find_vertex(fid)
in='';
while ~isequal(in,'Vertices') && ~feof(fid)
    in=fscanf(fid,'%s',1);
end;
function find_face(fid)
in='';
while ~isequal(in,'Triangles') && ~feof(fid)
    in=fscanf(fid,'%s',1);
end;