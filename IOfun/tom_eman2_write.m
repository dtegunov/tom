function tom_eman2_write(filename,data,flag)
%TOM_EMAN2_WRITE writes hdf5 format for eman2 and sparx
%
%   tom_eman2_write(filename,data)
%
%PARAMETERS
%
%  INPUT
%   filename            filename of the emdata      
%   data                matrix with or without header 
%                       which should be written   
%   flag                use '2dstack','volume' or 'volstack' (4d data stack of volumes)   
%                       in e2 a volume differs from a 2d stack
%  OUTPUT
%    
%
%EXAMPLE
%
% tom_eman2_write('out.hdf',tom_spheremask(ones(64,64,64),15),'volume');
%
%REFERENCES
%
%SEE ALSO
%   h5info,h5write,tom_eman2_read
%
%   created by FB 09/06/12
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

%remove preexisting file 2 allow overwrite 
if (exist(filename,'file'))
    unix(['rm ' filename]);
    pause(0.45);
end;

%location 
loc4img='/MDF/images/';

%get size of input data
sz=size(data);
if (length(sz)==2)
    sz=[sz 1];
end;


%write it!
if (strcmp(flag,'volume'))
    h5create(filename,[loc4img '0' '/image'],sz,'Datatype','single');
    h5write(filename,[loc4img  '0' '/image'],single(data));
    h5writeatt(filename,loc4img,'imageid_max',int32(0));
    [mea maxa mina stda]=tom_dev(data,'noinfo');
    wirte_eman_att([loc4img  '0/'],filename,[mea maxa mina stda sqrt(sum(sum(sum(data))))],sz,sz(3));
end;

if (strcmp(flag,'2dstack'))
    for i=1:size(data,3)
        h5create(filename,[loc4img num2str(i-1) '/image'],sz(1:2),'Datatype','single');
        h5write(filename,[loc4img  num2str(i-1) '/image'],single(data(:,:,i)));
        h5writeatt(filename,loc4img,'imageid_max',int32(size(data,3)-1));
        [mea maxa mina stda]=tom_dev(data,'noinfo');
        wirte_eman_att([loc4img  num2str(i-1)],filename,[mea maxa mina stda (sum(sum(sum(data)))).^2 ],sz,1);
    end;
end;


function wirte_eman_att(e2loca,filename,stat,sz,nz)


h5writeatt(filename,e2loca,'EMAN.DDD.imageid_max',int32(0));
h5writeatt(filename,e2loca,'EMAN.HostEndian','little');
h5writeatt(filename,e2loca,'EMAN.ImageEndian','little');
%h5writeatt(filename,e2loca,'EMAN.ImageEndian','big');
h5writeatt(filename,e2loca,'EMAN.apix_x',1.000000);
h5writeatt(filename,e2loca,'EMAN.apix_y',1.000000);
h5writeatt(filename,e2loca,'EMAN.apix_z',1.000000);
h5writeatt(filename,e2loca,'EMAN.changecount',int32(2));
h5writeatt(filename,e2loca,'EMAN.datatype',int32(7));
h5writeatt(filename,e2loca,'EMAN.is_complex',0);
h5writeatt(filename,e2loca,'EMAN.is_complex_ri',1);
h5writeatt(filename,e2loca,'EMAN.is_complex_x',0);
h5writeatt(filename,e2loca,'EMAN.maximum',stat(2));
h5writeatt(filename,e2loca,'EMAN.mean',stat(1));
h5writeatt(filename,e2loca,'EMAN.mean_nonzero',stat(1));
h5writeatt(filename,e2loca,'EMAN.minimum',stat(3));
h5writeatt(filename,e2loca,'EMAN.nx',int32(sz(1)));
h5writeatt(filename,e2loca,'EMAN.ny',int32(sz(2)));
h5writeatt(filename,e2loca,'EMAN.nz',int32(nz));
h5writeatt(filename,e2loca,'EMAN.sigma',stat(4));
h5writeatt(filename,e2loca,'EMAN.sigma_nonzero',stat(4));
h5writeatt(filename,e2loca,'EMAN.source_n',0);
h5writeatt(filename,e2loca,'EMAN.source_path',filename);
h5writeatt(filename,e2loca,'EMAN.square_sum',stat(5));



