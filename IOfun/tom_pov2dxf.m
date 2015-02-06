function total_faces=tom_pov2dxf(filename,dxf_filename,reduce_patch,matlab_surface_file)

% TOM_POV2DXF converts a .pov file created by chimera into a .dxf file
% which can be used by Bryce. TOM_POV2DXF can handle large files and allows
% surface reduction. 20000-30000 is a good number of faces for Bryce.
%
%   tom_pov2dxf(filename,dxf_filename,reduce_patch,matlab_surface_file)
%
% PARAMETERS
%   INPUT
%   filename        name of the .pov file
%   dxf_filename    name of the .dxf file 
%   reduce_patch    reduces the number of triangles (0<reduce_patch<=1).
%                   1.0 is 100 percent (no reduction).
%                   A number > 1 defines the number of faces which will be created.
%   matlab_surface_file Matlab .mat filename with surface object (optional).
%
%   OUTPUT
%   total_faces     total number of faces created
%
% EXAMPLE
%   tom_pov2dxf('blob.pov','blob.dxf',.5,'surface.mat');
%   load surface.mat
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
%   tom_pov2dxf_gui, tom_stl2surface, tom_surface2stl, tom_surface2dxf,
%   reducepatch, tom_surface2HxSurface, tom_HxSurface2surface
%
%   last change
%   07/13/08 SN
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

wait = waitbar(10,'Converting POV file to DXF file ...');
set(wait,'Name',['Create temporary files.']); drawnow;
t0=clock;
actual_time=0;

if nargin<4 no_surface_file=0; else no_surface_file=1;end;

x_all=[];
y_all=[];
z_all=[];
f1_all=[];
f2_all=[];
f3_all=[];
%xn_all=[]; % normals are not recomputed after reduction. We put that in
%when Matlab shows that feature. Maybe 2011:-)).
%yn_all=[];
%zn_all=[];
vertices=[];
faces=[];
normals=[];
nr_meshes=0;
fid = fopen(filename);
if fid==-1
    disp('cannot open .pov file.'); close(wait); total_faces=0; return;
end;
tic;
fseek(fid,0,1);
full_size=ftell(fid);
fseek(fid,0,-1);

mesh_idx=1;
patch_nr_of_meshes=1000;
total_nr_of_meshes=0;

create_dxf_file(dxf_filename);
% update_dxf_header(dxf_filename,6540,100001);
total_vertices=0;
total_faces=0;
nr_red_vertices=0;
all_nr_vertices=0;

while ~feof(fid)

    while ~feof(fid)
        find_mesh2(fid);
        [nr nr_faces x y z f1 f2 f3 xn yn zn]=read_mesh2(fid);
        x_all=[x_all; x];
        y_all=[y_all; y];
        z_all=[z_all; z];
        f1_all=[f1_all; f1+mesh_idx];
        f2_all=[f2_all; f2+mesh_idx];
        f3_all=[f3_all; f3+mesh_idx];
        %        xn_all=[xn_all; xn];
        %        yn_all=[yn_all; yn];
        %        zn_all=[zn_all; zn];
        mesh_idx=mesh_idx+size(x,1);
        %        face_idx=face_idx+size(f1,1);
        nr_meshes=nr_meshes+1;
        all_nr_vertices=all_nr_vertices+nr;
        %        disp([num2str(x(1)) ' ' num2str(y(1))]);
        %   end;
        %    if mod(nr_meshes,500)==0
        %    disp('.');
        %    end;

        if nr_meshes==patch_nr_of_meshes; break;end;
    end;
    current_postion=ftell(fid);
    full_time=num2str(round((full_size./current_postion).*etime(clock,t0)));
    full_percent_time=num2str(round(current_postion./full_size.*100));
    set(wait,'Name',[' ' num2str(round(etime(clock,t0))) ' s of in total ' full_time ' s. ' full_percent_time ' % done.']);
    waitbar(current_postion./full_size,wait);
    g.vertices=[x_all y_all z_all];
    g.faces=[f1_all f2_all f3_all];
    %    g.normals=[xn_all yn_all zn_all];

    x_all=[]; y_all=[]; z_all=[]; f1_all=[]; f2_all=[]; f3_all=[];
    total_nr_of_meshes=total_nr_of_meshes+nr_meshes;
    nr_meshes=0;
    mesh_idx=1;
    p_h=figure;
    p=patch(g);
    %   set(p,'NormalMode','manual');
    %   set(p,'VertexNormals',[xn_all yn_all zn_all ]);
    %   xn_all=[]; yn_all=[]; zn_all=[];

    try
        reducepatch(p,reduce_patch);
        red_vertices=get(p,'Vertices');
        red_faces=get(p,'Faces');
        %        red_normals=get(p,'VertexNormals'); % normals are not recomputed
        %        after reduction
    catch
        red_vertices=g.vertices;
        red_faces=g.faces;
        %        red_normals=get(p,'VertexNormals');
    end;
    %    normals=[normals; red_normals];

    red_faces=red_faces+nr_red_vertices;
    nr_red_vertices=nr_red_vertices+size(red_vertices,1);
    close(p_h);
    clear p; 

    write_patch_to_dxf(dxf_filename,red_vertices,red_faces);
    total_vertices=size(red_vertices,1)+total_vertices;
    total_faces=size(red_faces,1)+total_faces;

    if no_surface_file
        vertices=[vertices; red_vertices];
        faces=[faces; red_faces];
    end;
end;
fclose(fid);
set(wait,'Name',['Write file and clean-up.']); drawnow;

if no_surface_file
    surface.vertices=vertices;
    surface.faces=faces;
    save(matlab_surface_file,'surface');
end;
%surface.normals=normals;
%create_dxf_header(dxf_filename,total_vertices,total_faces);
attach_vertices_faces_dxf_file(dxf_filename,total_vertices,total_faces);
%close_dxf_file(dxf_filename);
%close(gcf);
% delete the temporary files
name=sprintf('%s.vertices.dxf',dxf_filename);
delete(name);
name=sprintf('%s.faces.dxf',dxf_filename);
delete(name);
close(wait);
disp(['Total number of faces created: ' num2str(total_faces)])

function [nr nr_faces x y z f1 f2 f3 xn yn zn]=read_mesh2(fid)
%        disp('read_mesh2')
[s]=fscanf(fid,'%s%s%s',3);
[nr]=fscanf(fid,'%u',1);
x=zeros([nr 1],'single');
y=zeros([nr 1],'single');
z=zeros([nr 1],'single');
find_bra(fid);
%        [s]=fscanf(fid,'%s',1);
%        [s]=fscanf(fid,'%c',7); % seems to be [s]=fscanf(fid,'%c',6); for linux
%       disp(['reading the vertex_vectors, nr:' num2str(nr)]);
for i=1:nr
    [x(i)]=fscanf(fid,'%f',1);
    [s]=fscanf(fid,'%c',2);
    y(i)=fscanf(fid,'%f',1);
    [s]=fscanf(fid,'%c',2);
    z(i)=fscanf(fid,'%f',1);
    s=fscanf(fid,'%c',3);
end;
xn=zeros([nr 1],'single');
yn=zeros([nr 1],'single');
zn=zeros([nr 1],'single');
find_bra(fid);
%        [s]=fscanf(fid,'%s%s%s',3);
%        [s]=fscanf(fid,'%c',7); % seems to be [s]=fscanf(fid,'%c',6); for linux
%      disp(['reading the normal_vectors, nr:' num2str(nr)]);
for i=1:nr
    [xn(i)]=fscanf(fid,'%f',1);
    [s]=fscanf(fid,'%c',2);
    yn(i)=fscanf(fid,'%f',1);
    [s]=fscanf(fid,'%c',2);
    zn(i)=fscanf(fid,'%f',1);
    s=fscanf(fid,'%c',3);
end;
[s]=fscanf(fid,'%s%s',2);
[nr_faces]=fscanf(fid,'%u',1);
find_bra(fid);
%[s]=fscanf(fid,'%c',8); % seems to be [s]=fscanf(fid,'%c',7); for linux
%     disp(['reading the faces_indices, nr:' num2str(nr)]);
f1=zeros([nr_faces 1],'uint32');
f2=zeros([nr_faces 1],'uint32');
f3=zeros([nr_faces 1],'uint32');
for i=1:nr_faces
    f1(i)=fscanf(fid,'%u',1);
    [s]=fscanf(fid,'%c',2);
    f2(i)=fscanf(fid,'%u',1);
    [s]=fscanf(fid,'%c',2);
    f3(i)=fscanf(fid,'%u',1);
    s=fscanf(fid,'%c',3);
end;


function find_mesh2(fid)
in='';
while ~isequal(in,'mesh2') & ~feof(fid)
    in=fscanf(fid,'%s',1);
end;
function find_bra(fid)
in='';
while ~isequal(in,'<') & ~feof(fid)
    in=fscanf(fid,'%c',1);
end;
function find_cket(fid)
in='';
while ~isequal(in,'>') & ~feof(fid)
    in=fscanf(fid,'%c',1);
end;

%            fseek(fid,0,-5);
function create_dxf_header(dxf_filename,nr_vertices,nr_faces)
name_dxf=dxf_filename;
name=sprintf('%s.header.dxf',dxf_filename);
when=datestr(now);
Header=sprintf(['999\ncreated by Matlab and the TOM toolbox at: ' when '\n0\nSECTION\n2\nENTITIES\n0\n']);
fid=fopen(name,'w');
fprintf(fid,Header);
fprintf(fid,'POLYLINE\n8\n%s\n66\n1\n70\n64\n71\n%d\n72\n%d\n62\n144\n0\n',name_dxf,nr_vertices,nr_faces);
fclose(fid);
function create_dxf_file(dxf_filename)
name_dxf=dxf_filename;
fid=fopen(name_dxf,'w'); if fid==-1 disp('cannot write DXF file. Permissions???'); end;
fclose(fid);
name_vertices=sprintf('%s.vertices.dxf',dxf_filename);
fid=fopen(name_vertices,'w');
fclose(fid);
name_faces=sprintf('%s.faces.dxf',dxf_filename);
fid=fopen(name_faces,'w');
fclose(fid);

%                function update_dxf_header(dxf_filename,aantalVertex,aantalFaces)
%                name=sprintf('%s.dxf',dxf_filename);
%                fid=fopen(name,'r+');
%                fseek(fid,59,0);
%                fprintf(fid,'POLYLINE\n8\n%s\n66\n1\n70\n64\n71\n%d\n72\n%d\n62\n144\n0\n',name,aantalVertex, (2*aantalFaces));
%                fclose(fid);

function write_patch_to_dxf(dxf_filename,vertices,faces)
name=dxf_filename;
name_vertices=sprintf('%s.vertices.dxf',dxf_filename);
fid=fopen(name_vertices,'a+');
for i=1:size(vertices,1)
    fprintf(fid,'VERTEX\n8\n%s\n',name);
    fprintf(fid,'10\n%.4f\n20\n%.4f\n30\n%.4f\n',vertices(i,1),vertices(i,2),vertices(i,3));
    fprintf(fid,'70\n192\n0\n');
end;
fclose(fid);
name_faces=sprintf('%s.faces.dxf',dxf_filename);
name=dxf_filename;
fid=fopen(name_faces,'a+');
for i=1:size(faces,1)
    %header
    fprintf(fid,'VERTEX\n8\n%s\n10\n0\n20\n0\n30\n0\n70\n128\n',name);
    %connections
    fprintf(fid,'71\n%d\n72\n%d\n73\n-%d\n',faces(i,1),faces(i,2),faces(i,3));
    %end
    fprintf(fid,'0\n');
    %header
    fprintf(fid,'VERTEX\n8\n%s\n10\n0\n20\n0\n30\n0\n70\n128\n',name);
    %connections
    fprintf(fid,'71\n%d\n72\n%d\n73\n-%d\n',faces(i,3),faces(i,2),faces(i,1));
    %end
    fprintf(fid,'0\n');
end;
fclose(fid);

function close_dxf_file(dxf_filename)
name=dxf_filename;
fid=fopen(name,'a+');
fprintf(fid,'SEQEND\n8\n%s\n0\n',name);
fprintf(fid,'ENDSEC\n0\nEOF\n');
fclose(fid);

function attach_vertices_faces_dxf_file(dxf_filename,nr_vertices,nr_faces)
name_dxf=dxf_filename;
name_vertices=sprintf('%s.vertices.dxf',dxf_filename);
fid=fopen(name_vertices,'r+');
fid_dxf=fopen(name_dxf,'w+');
when=datestr(now);
Header=sprintf(['999\ncreated by Matlab and the TOM toolbox at: ' when '\n0\nSECTION\n2\nENTITIES\n0\n']);
fprintf(fid_dxf,Header);
fprintf(fid_dxf,'POLYLINE\n8\n%s\n66\n1\n70\n64\n71\n%d\n72\n%d\n62\n144\n0\n',name_dxf,nr_vertices,2.*nr_faces);

while feof(fid)==0
    [s]=fscanf(fid,'%s',1);
    if ~isequal(s,'')
        fprintf(fid_dxf,'%s\n',s);
    end;
end;
fclose(fid);
name_faces=sprintf('%s.faces.dxf',dxf_filename);
fid=fopen(name_faces,'r+');
while feof(fid)==0
    [s]=fscanf(fid,'%s',1);
    if ~isequal(s,'')
        fprintf(fid_dxf,'%s\n',s);
    end;
end;

fprintf(fid_dxf,'SEQEND\n8\n%s\n0\n',name_dxf);
fprintf(fid_dxf,'ENDSEC\n0\nEOF\n');

fclose(fid);
fclose(fid_dxf);

% Definition POLYLINE, 2002
% Polyline group codes
% Group codes       Description
% 100               Subclass marker (AcDb2dPolyline or AcDb3dPolyline)
% 66                Obsolete; formerly an “entities follow flag” (optional; ignore if present)
% 10                DXF: always 0
%                   APP: a “dummy” point; the X and Y values are always 0, and the Z value
%                   is the polyline’s elevation (in OCS when 2D, WCS when 3D)
% 20                DXF: always 0
% 30                DXF: polyline’s elevation (in OCS when 2D, WCS when 3D)
% 39                Thickness (optional; default = 0)
% 70                Polyline flag (bit-coded); default is 0:
%                   1 = This is a closed polyline (or a polygon mesh closed in the M
%                   direction)
%                   2 = Curve-fit vertices have been added
%                   4 = Spline-fit vertices have been added
%                   8 = This is a 3D polyline
%                   16 = This is a 3D polygon mesh
%                   32 = The polygon mesh is closed in the N direction
%                   64 = The polyline is a polyface mesh
%                   128 = The linetype pattern is generated continuously around the
%                   vertices of this polyline
% 40                Default start width (optional; default = 0)
% 41                Default end width (optional; default = 0)
% 71                Polygon mesh M vertex count (optional; default = 0)
% 72                Polygon mesh N vertex count (optional; default = 0)
% 73                Smooth surface M density (optional; default = 0)
% 74                Smooth surface N density (optional; default = 0)

%Group codes        Description
%100                Subclass marker (AcDbVertex)
%100                Subclass marker (AcDb2dVertex or AcDb3dPolylineVertex)
%10                 Location point (in OCS when 2D, and WCS when 3D)
%                   DXF: X value; APP: 3D point
%20, 30             DXF: Y and Z values of location point (in OCS when 2D, and WCS when
%                   3D)
%40                 Starting width (optional; default is 0)


% Definition 3DFACE, 2002
% Group codes   Description
% 100           Subclass marker (AcDbFace)
% 10            First corner (in WCS)
%               DXF: X value; APP: 3D point
% 20, 30        DXF: Y and Z values of first corner (in WCS)
% 11            Second corner (in WCS)
%               DXF: X value; APP: 3D point
% 21, 31        DXF: Y and Z values of second corner (in WCS)
% 12            Third corner (in WCS)
%               DXF: X value; APP: 3D point
% 22, 32        DXF: Y and Z values of third corner (in WCS)
% 13            Fourth corner (in WCS). If only three corners are entered, this is the
%               same as the third corner
%               DXF: X value; APP: 3D point
% 23, 33        DXF: Y and Z values of fourth corner (in WCS)
% 70            Invisible edge flags (optional; default = 0):
%               1 = First edge is invisible
%               2 = Second edge is invisible
%               4 = Third edge is invisible
%               8 = Fourth edge is invisible
