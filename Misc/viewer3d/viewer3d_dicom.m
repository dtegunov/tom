function [voxelvolume,scales]=viewer3d_dicom(dirname)
% This function viewer3d_dicom reads a folder containing the 
% dicom slice files of one volume and returns it as a matlab volume.
%
% Examples:
%
% 1. I = viewer3d_dicom;
% 2. I = viewer3d_dicom('c:\patient\dicom');
% 3. [I,scales] = viewer3d_dicom;
%
% (I contains the imagevolume, scales the pixel dimensions)

% Check if function is called with folder name
if(exist('dirname','var')==0)
    [filename, dirname] = uigetfile( {'*.dcm;*.dicom', 'Dicom Files'; '*.*', 'All Files (*.*)'}, 'Select a dicom file');
    fullfn=[dirname filename];
    voxelvolume=squeeze(dicomread(fullfn));
    if(ndims(voxelvolume)>2)
        scales=[1 1 1]; 
        warning('viewer3d_dicom:file','voxel scales are not known');
        return;
    end
end
  

% Read directory for Dicom File Series
filelist = dir([dirname '\*.*']);
dicomfilelist=getdicomfilelist(filelist,dirname);
if(isempty(dicomfilelist)), 
    error('viewer3d_dicom:files','No dicom files found');
end

% Convert Dicom images to Voxel Volume
sizeh=length(dicomfilelist);
for i=1:sizeh;
    filename=filelist(dicomfilelist(i,2)).name;
    link=[dirname '\' filename];
    info = dicominfo(link);
    dicomimage = dicomread(info);
    if(i==1)
        [sizel,sizew]=size(dicomimage);
        voxelvolume=zeros([sizel sizew sizeh],class(dicomimage));
        SlicezA=info.ImagePositionPatient;
        SliceSpacing=info.SliceThickness;
    end
    if(i==2)
        SlicezB=info.ImagePositionPatient;
        SliceSpacing=abs(SlicezB(3)-SlicezA(3));
    end
    voxelvolume(:,:,i)=dicomimage;
end

% Pixel dimensions
scales=[info.PixelSpacing;SliceSpacing]';

function dicomfilelist=getdicomfilelist(filelist,dirname)
dicomfilelist=[]; m=0;
for i=1:length(filelist)
    if(file_is_dicom([dirname '\' filelist(i).name]))
        % Get the file number from the filename
        filenumstr=filelist(i).name; filenumstr(uint8(filenumstr)<48)=[]; filenumstr(uint8(filenumstr)>57)=[];
        m=m+1;
        dicomfilelist(m,1:2)=[str2num(filenumstr) i];
    end
end
dicomfilelist=sortrows(dicomfilelist);

function isdicom=file_is_dicom(filename)
isdicom=false;
try
    fid = fopen(filename, 'r');
    status=fseek(fid,128,-1);
    if(status==0)
        tag = fread(fid, 4, 'uint8=>char')';
        isdicom=strcmpi(tag,'DICM');
    end
    fclose(fid);
catch me
end

