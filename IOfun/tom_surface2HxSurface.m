function tom_surface2HxSurface(matlab_surface,hx_filename)

% TOM_SURFACE2HXSURFACE converts a matlab surface structure to a HxSurface file
% in ASCII format which can be used by Amira.
%
%   tom_surface2HxSurface(matlab_surface,hx_filename)
%
% PARAMETERS
%   INPUT
%   matlab_surface surface-structure
%   hx_filename    name of the HxSurface file
%
%   OUTPUT
%   Hx-file        file in HxSurface format
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
%   %And write it to an HxSurface-file:
%   tom_surface2HxSurface(surface,'simplex.surf')
%
% REFERENCES
%
% SEE ALSO
%   tom_pov2dxf, tom_pov2dxf_gui, tom_surface2stl, tom_surface2dxf, patch,
%   reducepatch, tom_HxSurface2surface
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
create_hx_file_ascii(matlab_surface,hx_filename)
write_hx_file_ascii(matlab_surface,hx_filename);

function write_hx_file_ascii(matlab_surface,hx_filename)

fid=fopen(hx_filename,'a'); if fid==-1 error('cannot write Hx file. Permissions???'); end;

Header=sprintf(['Vertices ' num2str(size(matlab_surface.vertices,1)) '\n']);
fprintf(fid,Header);

for nr_vertex=1:size(matlab_surface.vertices,1)
    v=matlab_surface.vertices(nr_vertex,:);
    vertex=sprintf(' %f %f %f \n',v(1), v(2) ,v(3));
    fprintf(fid,vertex);
end;
loop=sprintf(['Patches 1' '\n{\n']);
fprintf(fid,loop);
Header=sprintf(['Triangles ' num2str(size(matlab_surface.faces,1)) '\n']);
fprintf(fid,Header);
for nr_faces=1:size(matlab_surface.faces,1)
    v=matlab_surface.faces(nr_faces,:);
    face=sprintf(' %i %i %i \n',v(1), v(2) ,v(3));
    fprintf(fid,face);
end;

loop=sprintf(['}' '\n']); %close the file
fprintf(fid,loop);
fclose(fid);

function create_hx_file_ascii(matlab_surface,hx_filename)
fid=fopen(hx_filename,'w'); if fid==-1 error('cannot write Hx file. Permissions???'); end;
Header=sprintf(['# HyperSurface 0.1 ASCII\n']);
fprintf(fid,Header);
when=datestr(now);
Header=sprintf(['# Created by Matlab and the TOM toolbox at: ' when '\n']);
fprintf(fid,Header);
Header=sprintf(['Parameters {\n']);
fprintf(fid,Header);
Header=sprintf([' Filename "' hx_filename '"\n']);
fprintf(fid,Header);
Header=sprintf(['}\n']);
fprintf(fid,Header);
fclose(fid);