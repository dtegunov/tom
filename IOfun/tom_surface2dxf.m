function tom_surface2dxf(matlab_surface,dxf_filename)

% TOM_SURFACE2DXF converts a matlab surface structure to a .dxf file
% which can be used by Bryce.
%
%   tom_surface2dxf(matlab_surface,dxf_filename)
%
% PARAMETERS
%   INPUT
%   matlab_surface  surface-structure
%   dxf_filename    name of the .dxf file 
%
%   OUTPUT
%   dxf-file        file in DXF format
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
%   %And write it to an .dxf-file:
%   tom_surface2dxf(surface,'simplex.dxf')
%
%REFERENCES
%
% SEE ALSO
%   tom_pov2dxf, tom_pov2dxf_gui, tom_surface2stl, tom_stl2surface, patch,
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

nr_vertices=size(matlab_surface.vertices,1);
nr_faces=size(matlab_surface.faces,1);
create_dxf_header(dxf_filename,nr_vertices,nr_faces);
write_patch_to_dxf(dxf_filename,matlab_surface.vertices,matlab_surface.faces);
close_dxf_file(dxf_filename)

%end;
function create_dxf_header(dxf_filename,nr_vertices,nr_faces)
when=datestr(now);
Header=sprintf(['999\ncreated by Matlab and the TOM toolbox at: ' when '\n0\nSECTION\n2\nENTITIES\n0\n']);
fid=fopen(dxf_filename,'w');
fprintf(fid,Header);
fprintf(fid,'POLYLINE\n8\n%s\n66\n1\n70\n64\n71\n%d\n72\n%d\n62\n144\n0\n',dxf_filename,nr_vertices,nr_faces);
fclose(fid);

function write_patch_to_dxf(dxf_filename,vertices,faces)
fid=fopen(dxf_filename,'a+');
for i=1:size(vertices,1)
    fprintf(fid,'VERTEX\n8\n%s\n',dxf_filename);
    fprintf(fid,'10\n%.4f\n20\n%.4f\n30\n%.4f\n',vertices(i,1),vertices(i,2),vertices(i,3));
    fprintf(fid,'70\n192\n0\n');
end;
for i=1:size(faces,1)
    %header
    fprintf(fid,'VERTEX\n8\n%s\n10\n0\n20\n0\n30\n0\n70\n128\n',dxf_filename);
    %connections
    fprintf(fid,'71\n%d\n72\n%d\n73\n-%d\n',faces(i,1),faces(i,2),faces(i,3));
    %end
    fprintf(fid,'0\n');
    %header
    fprintf(fid,'VERTEX\n8\n%s\n10\n0\n20\n0\n30\n0\n70\n128\n',dxf_filename);
    %connections
    fprintf(fid,'71\n%d\n72\n%d\n73\n-%d\n',faces(i,3),faces(i,2),faces(i,1));
    %end
    fprintf(fid,'0\n');
end;
fclose(fid);

function close_dxf_file(dxf_filename)
fid=fopen(dxf_filename,'a+');
fprintf(fid,'SEQEND\n8\n%s\n0\n',dxf_filename);
fprintf(fid,'ENDSEC\n0\nEOF\n');
fclose(fid);
