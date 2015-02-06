function stack=tom_eman2_read(filename)
%TOM_EMAN2_READ reads hdf5 format from eman2 and sparx
%
%   stack=tom_eman2_read(filename)
%
%PARAMETERS
%
%  INPUT
%   filename            filename of the emdata      
%  
%  
%  OUTPUT
%    data               data in memory  
%
%EXAMPLE
%
%vol = tom_eman2_read('start.hdf');
%
%REFERENCES
%
%SEE ALSO
%   h5info,h5read,tom_spiderrread
%
%   created by FB 09/05/12
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

%default locations
loc_num_of_entry='/MDF/images/';
att_name_num_of_entry='imageid_max';
loc_of_img='/image';

%get num_of_entries
num_of_entries = h5readatt(filename,loc_num_of_entry,att_name_num_of_entry)+1;

%get data dimension
tmp_img=h5read(filename,[loc_num_of_entry '0' loc_of_img]);
dim_of_entry=size(tmp_img);

%allocate memory for reading
if (num_of_entries > 1)
    stack=zeros([dim_of_entry num_of_entries],'single');
else
    stack=zeros(dim_of_entry,'single');
end;

for i=1:num_of_entries
    if (length(dim_of_entry)>2)
        stack(:,:,:,i)=h5read(filename,[loc_num_of_entry num2str(i-1) loc_of_img]);
    else
        stack(:,:,i)=h5read(filename,[loc_num_of_entry num2str(i-1) loc_of_img]);
    end;
 end;

stack=tom_emheader(stack);



