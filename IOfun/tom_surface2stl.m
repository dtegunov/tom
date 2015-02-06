function tom_surface2stl(matlab_surface,stl_filename,format,matlab_normals)

% TOM_SURFACE2STL converts a matlab surface structure to a .stl file
% which can be used by Amira.
%
%   tom_surface2stl(matlab_surface,stl_filename)
%
% PARAMETERS
%   INPUT
%   matlab_surface  surface-structure
%   stl_filename    name of the .stl file 
%   format          file format: 'ascii' (default) or 'binary'
%   matlab_normals  given normals, a vector with the size of faces
%
%   OUTPUT
%   stl-file        file in STL format
%
% An STL (“StereoLithography”) file is a triangular representation
% of a 3-dimensional surface geometry. The surface is tessellated
% or broken down logically into a series of small triangles (facets).
% Each facet is described by a perpendicular direction and three
% points representing the vertices (corners) of the triangle. 
%
% STL-ASCII-Format:
% solid name
%   facet normal -0.577349 -0.577349 -0.577354
%       outer loop
%           vertex   137.660370 251.773575 286.188660
%           vertex   135.849045 253.584900 286.188660
%           vertex   137.660370 253.584900 284.377350
%       endloop
%   endfacet
% endsolid name
% 
% STL-Binary-Format:
% Bytes     Data type               Description
%   80      ASCII                   Header No information needed
%   4       unsigned long integer   Number of facets
%   4       float i                 normal-vector
%   4       float j                 normal-vector
%   4       float k                 normal-vector
%   4       float x                 vertex_1-X
%   4       float y                 vertex_1-Y
%   4       float z                 vertex_1-Z
%   4       float x                 vertex_2-X
%   4       float y                 vertex_2-Y
%   4       float z                 vertex_2-Z
%   4       float x                 vertex_3-X
%   4       float y                 vertex_3-Y
%   4       float z                 vertex_3-Z
%   2       unsigned integer        not used, default 0
%
% EXAMPLE
%   %Let's define a simplex surface:
%   surface.vertices=...
%   [0 0 0; ...
%    1 0 0; ...
%    1 1 0; ...
%    1 1 1 ];
%   surface.faces=...
%   [1 2 3; ...
%    1 2 4; ...
%    1 3 4; ...
%    2 3 4 ];
%   %Show me the surface:
%   p=patch('Faces',surface.faces,'Vertices',surface.vertices,'FaceVertexCData'...
%   ,ones(size(surface.faces,1),3).*.8,'FaceColor','flat');
%   view(3);grid on; axis equal;
%   %And write it to an .stl-file:
%   tom_surface2stl(surface,'simplex.stl')
%
% REFERENCES
%
% SEE ALSO
%   tom_pov2dxf, tom_pov2dxf_gui, tom_surface2stl, tom_surface2dxf, patch,
%   reducepatch
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

if (~(isfield(matlab_surface,'vertices') && isfield(matlab_surface,'faces')))
    error( 'matlab_surface must have a structure format' );
end

if nargin<3
    format='ascii';
end;

if isequal(format,'ascii')
    create_stl_file_ascii(matlab_surface,stl_filename)
    if exist('matlab_normals','var') 
        write_stl_file_ascii(matlab_surface,stl_filename,matlab_normals);
    else
        write_stl_file_ascii(matlab_surface,stl_filename);
    end;
else
    create_stl_file_binary(matlab_surface,stl_filename);
    if exist('matlab_normals','var')     
        write_stl_file_binary(matlab_surface,stl_filename,matlab_normals);
    else
        write_stl_file_binary(matlab_surface,stl_filename);
    end;        
end;

function write_stl_file_ascii(matlab_surface,stl_filename,matlab_normals);

fid=fopen(stl_filename,'a'); if fid==-1 error('cannot write STL file. Permissions???'); end;

for nr_faces=1:size(matlab_surface.faces,1)
    faces=matlab_surface.faces(nr_faces,:);
    v_1=matlab_surface.vertices(faces(1),:);
    v_2=matlab_surface.vertices(faces(2),:);
    v_3=matlab_surface.vertices(faces(3),:);
    if exist('matlab_normals','var') 
        normal=matlab_normals(nr_faces,:);
    else
        normal=calc_normal(v_1,v_2,v_3);
    end;
    normal=sprintf('facet normal %f %f %f \n', normal(1), normal(2), normal(3));
    fprintf(fid,normal);
    loop=sprintf(['outer loop' '\n']);
    fprintf(fid,loop);
    for nr_vertex=1:3;
        v=matlab_surface.vertices(faces(nr_vertex),:);
        vertex=sprintf(' vertex %f %f %f \n',v(1) ,v(2) ,v(3));
        fprintf(fid,vertex);
    end;
    loop=sprintf([' endloop' '\n' 'endfacet' '\n']);
    fprintf(fid,loop);
end;
loop=sprintf(['endsolid' '\n']);
fprintf(fid,loop);
fclose(fid);

function write_stl_file_binary(matlab_surface,stl_filename,matlab_normals);

fid=fopen(stl_filename,'ab+'); if fid==-1 error('cannot write STL file. Permissions???'); end;

for nr_faces=1:size(matlab_surface.faces,1)
    faces=matlab_surface.faces(nr_faces,:);
    v_1=matlab_surface.vertices(faces(1),:);
    v_2=matlab_surface.vertices(faces(2),:);
    v_3=matlab_surface.vertices(faces(3),:);
    if exist('matlab_normals','var') 
        normal=matlab_normals(nr_faces,:);
    else
        normal=calc_normal(v_1,v_2,v_3);
    end;
    fwrite(fid,normal,'float');
    for nr_vertex=1:3;
        v=matlab_surface.vertices(faces(nr_vertex),:);
        fwrite(fid,v,'float');
    end;
    fwrite(fid,0,'short');
end;
fclose(fid);

function create_stl_file_ascii(matlab_surface,stl_filename)
fid=fopen(stl_filename,'w'); if fid==-1 error('cannot write STL file. Permissions???'); end;
Header=sprintf(['solid ' stl_filename ' (' num2str(size(matlab_surface.faces,1)) ' triangles, ' num2str(size(matlab_surface.vertices,1)) ' nodes) \n']);
fprintf(fid,Header);
fclose(fid);

function create_stl_file_binary(matlab_surface,stl_filename)
fid=fopen(stl_filename,'wb'); if fid==-1 error('cannot write STL file. Permissions???'); end;
fwrite(fid,zeros(80,1),'uchar');
fwrite(fid,size(matlab_surface.faces,1),'unsigned long');
fclose(fid);

function normal=calc_normal(v_1,v_2,v_3)
d_1=v_2-v_1;
d_2=v_3-v_1;
d_3=cross(d_1,d_2);
normal=d_3./sqrt(sum(d_3.*d_3));