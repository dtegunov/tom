function tom_k2_sum2stack(pickListSum,docFileSum,frameFolders,fpickListFrames,fdocFileFrames,verbose)
%tom_k2_sum2stack transforms a sumed stack back 2 single frames
%
%    tom_k2_sum2stack(pickListSum,docFileSum,frameFolders,fpickListFrames,fdocFileFrames)
%
%PARAMETERS
%
%  INPUT
%   pickListSum            pickList with corrdinates on sumed stack micrographs
%   docFileSum             ('') xmipp-doc file with the angles from sumed stack
%   frameFolders           folder where the single frames (from Sumed stack) are located
%   fpickListFrames        filename of picklist for the frames
%   pickListPerFrame       (1) 1 picklist 4 every frame  
%   fdocFileFrames         ('') fllename of the docfile for the frames 
%   docFilePerFrame        (1)  1 doc file 4 every frame
%   verbose                (1) verbose flag    
%
%EXAMPLE
%   
% basep='/fs/pool/pool-titan5/K220S/data/0001_150313/'
% for i=1:10
%     frames4sp=[1 2 3]+(3*(i-1));
%     frameFolders{1}=[basep 'lowsp' num2str(frames4sp(1)) 'to'  num2str(frames4sp(end))];
% end;
% 
%  tom_k2_sum2stack(pickListSum,docFileSum,frameFolders,fpickListFrames,fdocFileFrames,verbose)
%
%NOTE
%
% 1 pickListSum and docFileSum have 2 have the same lenght and have to
%   correspond by index
% 2 Number of frames has to be consistent over the dataset
%
%SEE ALSO
%   tom_k2_average_stack
%
%   created by FB 03/18/13
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



if (nargin < 5)
    fdocFileFrames='';
end;



if (nargin < 6)
    verbose=1;
end;


load(pickListSum);
OrgAlign=align2d;

OrgDoc=tom_xmippdocread(docFileSum);

if (size(align2d,2)~=length(OrgDoc))
    error('docfile and picklist differ in length');
end;


%generate new align2d 4 frames
[pathOutPl,nameOutPl,extOutPl]=fileparts(fpickListFrames);

% new_alg=[];
% for i=1:length(frameFolders)
%     new_alg=cat(2,new_alg,align2d);
% end;

zz=1;
for ii=1:length(frameFolders)
    align2d=OrgAlign;
    [~,frame_name]=fileparts(frameFolders{ii});
    for i=1:size(OrgAlign,2)
        [tmp_path tmp_filename tmp_ext]=fileparts(OrgAlign(1,i).filename);
        align2d(1,i).filename=[frameFolders{ii} '_corr/' tmp_filename tmp_ext];
        %new_alg(1,zz)=align2d(1,i);
        zz=zz+1;
    end;
    save([pathOutPl '/' nameOutPl '_' frame_name  extOutPl],'align2d');
end;
%save(fpickListFrames,'new_alg');
disp('picklists done!');


for ii=1:length(frameFolders)
    doc=OrgDoc;
    [~,frame_name]=fileparts(frameFolders{ii});
    for i=1:length(OrgDoc)
        [tmp_path tmp_filename tmp_ext]=fileparts(OrgDoc(i).name);
        tmp_path=fileparts(tmp_path);
        [pref,num]=strtok(tmp_filename,'_');
        num=strrep(num,'_low_','');
        doc(i).name=[tmp_path '/frames/' frame_name '/parts/' pref '_' frame_name '_' num tmp_ext];
        %new_alg(1,zz)=align2d(1,i);
       
    end;
    warning off; mkdir([tmp_path '/frames/' frame_name]); warning on;
    tom_xmippdocwrite([tmp_path '/frames/' frame_name '/ang_doc.doc'],doc);
    unix(['cp ' pathOutPl '/' nameOutPl '_' frame_name  extOutPl ' ' tmp_path '/frames/' frame_name]);
end;

for ii=1:length(frameFolders)
    doc=OrgDoc;
    [~,frame_name]=fileparts(frameFolders{ii});
    [tmp_path,tmp_filename,tmp_ext]=fileparts(OrgDoc(1).name);
    [pref,num]=strtok(tmp_filename,'_');
    tmp_path=fileparts(tmp_path);
    NamePickList=[tmp_path '/frames/' frame_name '/' nameOutPl  '_' frame_name '.mat'];  
    basePath=[tmp_path '/frames/' frame_name '/parts/' pref '_' frame_name '_' ];
    SelName=[tmp_path '/frames/' frame_name '/parts.sel'];
    rad=160;
    tom_av2_xmipp_picklist2stack(NamePickList,basePath,SelName,rad,'gradient&mean0+1std-bg',1,[128 128],1);
    disp(' ');
    
end;


disp(' ');

function disp_it(message,verbose,log)

if (nargin < 3)
    log='';
end;

if (verbose==1)
    disp(message);
end;

if (isempty(log)==0)
    fid_log=fopen(log,'a');
    fprintf(fid_log,'%s\n',message);
    fclose(fid_log);
end;


