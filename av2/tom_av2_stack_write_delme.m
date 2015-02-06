function sucess_vector=tom_av2_stack_write(stack,outfilestruct,positions,save_outfile_flag)
%TOM_AV2_STACK_WRITE writes a stack in variable formats to Harddisk
%
%   tom_av2_stack_write(stack,outfilestruct,positions)
%
%PARAMETERS
%
%  INPUT
%   stack               stack as 3d matrix
%   outalign            output align struct
%   outfilestruct       structure for fileoutput
%                       outfilestruct.method             file structure '1stack' 'Nstacks' 'singleFiles'  
%                       outfilestruct.num_of_entries     number of entries
%                       outfilestruct.path               path for output             
%                       outfilestruct.folder_num_offset  folder numbering offset use -1 for empty 
%                       outfilestruct.part_nr_offset     (0) particle numbering offset 
%                       outfilestruct.filename           filename 
%                       outfilestruct.ext                extension     
%                       outfilestruct.fileformat         used fileformat (em,spider,imagic)
%                                               
%
%
%   positions            positions the particles in the stack should be
%                        written
%                        position can ether be a starting number or
%                        vector which has the same lenth as the stack
%                        or a flag 'append' (slow).
%                          
%  folder_num            folder numbering sceme 0 1 (default)  5 
%                        use -1 to start without number
%                                                
%                        
%                        
%
%
%  OUTPUT
%   error_pos           positions with writing error
%
%
%
%EXAMPLE
% 
%stack=rand(80,80,50);
%
%outfilestruct.method='singleFiles'; 
%outfilestruct.num_of_entries=10; %just for the show ...use 10000 for real application hurz!
%outfilestruct.path='part_fold';
%outfilestruct.folder_num_offset=1;
%
%outfilestruct.filename='part_';
%outfilestruct.ext='.em';
%outfilestruct.fileformat='em';
%outfilestruct.size=[80 80 50];
%
%tom_av2_stack_write(stack,outfilestruct,1);
%
% %append the stack on the harddisk
% tom_av2_stack_write(stack,outfilestruct,'append');
%
%
% write the stack on user defined positions ... useful for overwriting or
% or for preserving the correct order
% positions=200:249;
% tom_av2_stack_write(stack,outfilestruct,positions);
%
%   
%
%REFERENCES
%
%SEE ALSO
%   TOM_AV2_STACKBROWSER
%
%   created by FB(eckster) 06/16/08
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


if nargin < 2
    outfilestruct.method='singleFiles'; 
    outfilestruct.num_of_entries=10000;
    outfilestruct.path='part_fold';
    outfilestruct.filename='part_';
    outfilestruct.ext='.em';
    outfilestruct.fileformat='em';
    outfilestruct.part_nr_offset=0;
end;


if nargin < 3
    positions='append';
end;

if nargin < 4
    save_outfile_flag=0;
end;



if strcmp(outfilestruct.method,'singleFiles')==1

    %transfer variables
    folder_basename=outfilestruct.path;
    file_name=outfilestruct.filename;
    file_extension=outfilestruct.ext;
    files_per_folder=outfilestruct.num_of_entries;
    fileformat=outfilestruct.fileformat;
    folder_num_offset=outfilestruct.folder_num_offset;

    if ( not(max(size(positions))==1 || max(size(positions))==size(stack,3))  && strcmp(positions,'append')==0 )
        error('inproper writing position !!');
    end;

    if (strcmp(positions,'append')==1)
        positions=find_append_position(outfilestruct);
    end;

    if (max(size(positions))==1)
        positions=positions(1):positions(1)+size(stack,3);
    end;

    if (folder_basename(end)=='/')
       folder_basename=folder_basename(1:end-1); 
    end;

    sucess_vector=ones(size(stack,3),1);
    
    
    for i=1:size(stack,3)
        
        folder_num=floor((positions(i)-1)./files_per_folder)+folder_num_offset;
        
        if (folder_num==-1)
            folder_num='';
        end;
        
        if (mod(positions(i),files_per_folder) == 1 || i==1 || exist([folder_basename num2str(folder_num)],'dir')==0 )
            warning off; mkdir([folder_basename num2str(folder_num)]); warning on;
        end;
        
        name_string=[folder_basename num2str(folder_num)  '/' file_name num2str(positions(i))  file_extension ];
       
        
        sucess=save_write(name_string,stack(:,:,i),fileformat);
       
        if (sucess==0)
            %warning(['error writing ' name_string]);
            pause(1);
            sucess(i)=positions(i);
        end;
        
        
    end;

   
    if (save_outfile_flag==1)
        if (length(outfilestruct.size)==2)
            outfilestruct.size=[outfilestruct.size 1];
        end;
        save([folder_basename '.mat'],'outfilestruct');
    end;
    
    return;
    
end;


% if strcmp(outfilestruct.method,'singleFiles')==1
%     
% end;








function positions=find_append_position(outfilestruct)

[a b c]=fileparts(outfilestruct.path);
if isempty(a)
    p_path='./';
else
    p_path=a;
end;
dd=dir([p_path '/' b '*']);

if (isempty(dd)==0)
    for i=1:size(dd,1)
        in_cell{i}=dd(i).name;
    end;
    max_num=tom_get_max_numeric(in_cell);

    dd=dir([p_path '/' b  num2str(max_num) '/' '*' outfilestruct.ext]);
    if (isempty(dd)==0)
        for i=1:size(dd,1)
            in_cell{i}=dd(i).name;
        end;
        positions=tom_get_max_numeric(in_cell);
    else
        positions=1;
    end;
else
    positions=1;
end;


function sucess=save_write(name_string,out,fileformat)

sucess=0;
for i=1:10
    try
        switch fileformat
            case 'em'
                tom_emwrite(name_string,out);
            case 'spider'
                tom_spiderwrite(name_string,out);
            otherwise
                error('Fileformat not implemented !!!');
        end;
        sucess=1;
        return;
    catch
        pause(3);
    end;
end;




