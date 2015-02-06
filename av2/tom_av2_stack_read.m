function [stack positions_out]=tom_av2_stack_read(outfilestruct,positions,flag)
%TOM_AV2_STACK_READ reads a stack from Harddisk
%
%   [stack positions]=tom_av2_stack_read(outfilestruct,positions,flag)
%
%PARAMETERS
%
%  INPUT
%   outfilestruct       structure for fileoutput
%                       outfilestruct.method             file structure '1stack' 'Nstacks' 'singleFiles'  
%                       outfilestruct.num_of_entries     number of entries
%                       outfilestruct.path               path for output             
%                       outfilestruct.filename           filename 
%                       outfilestruct.ext                extension     
%                       outfilestruct.fileformat         used fileformat (em,spider,imagic)
%                                               
%
%
%   positions            positions the particles in the stack should be
%                        read
%                        position can ether be a starting number or
%                        vector which has the same lenth as the stack
%                        or a flag 'all' 
%
%   flag                 interpretation of the positions vector
%                           'start_and_length' vector with 2 entries goes through the entries until stack has reached required length  
%                           'explicit' [1 2 3 4 5] vector with all particle numbers the particles with the exact number 
%                           'start_and_length_by_position' vector with 2 entries goes through the entries until particle number reached required length    
%                                                           ==> stack can be shorter.    
%
%
%  OUTPUT
%   stack         
%
%EXAMPLE
% 
%
%
%outfilestruct.method='singleFiles'; 
%outfilestruct.num_of_entries=10; %just for the show ...use 10000 for real application hurz!
%outfilestruct.path='part_fold';
%outfilestruct.filename='part_';
%outfilestruct.ext='.em';
%outfilestruct.fileformat='em';
%
%stack=rand(80,80,50);
%tom_av2_stack_write(stack,outfilestruct,1);
%
%stack=tom_av2_stack_read(outfilestruct,[5 30],'start_and_length');
%
% %append the stack on the harddisk
% 
%
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

if nargin < 3
    flag='explicit';
end;
    

if strcmp(outfilestruct.method,'singleFiles')==1

    filesperfolder=outfilestruct.num_of_entries;
    folderbasename=outfilestruct.path;
    filename=outfilestruct.filename;
    file_extension=outfilestruct.ext;
    fold_num=floor((positions(1)-1)./filesperfolder)+1;
    fileformat=outfilestruct.fileformat;

    if (folderbasename(end)=='/')
        folderbasename=folderbasename(1:end-1);
    end;

    st=tom_reademheader([folderbasename num2str(fold_num) '/' filename num2str(positions(1)) '.em']);

    switch flag
        case 'start_and_length'
            req_length=positions(2);
            positions=positions(1):(positions(1)+20*positions(2));
        case 'explicit'
            req_length=length(positions);
        case 'start_and_length_by_position'
            req_length=positions(2);
            positions=positions(1):(positions(1)+positions(2));
        otherwise
            error('flag not implemented !!!');
    end;


    stack=zeros(st.Header.Size(1),st.Header.Size(2),req_length);

    zz=1;
    for i=1:length(positions)
        fold_num=floor((positions(i)-1)./filesperfolder)+1;
        name_str=[folderbasename num2str(fold_num) '/' filename num2str(positions(i)) file_extension];
        try
            switch fileformat
                case 'em'
                    tmp=tom_emreadc(name_str);
                case 'spider'
                    tmp=tom_siderread(name_str);
                otherwise
                    error('Fileformat not implemented !!!');
            end;
        catch
            if strcmp(flag,'explicit')
                warning([name_str ' not readable']);
            end;
            continue;
        end;
        positions_out(zz)=positions(i);
        stack(:,:,zz)=tmp.Value;
       

        if (zz==req_length)
            if (zz < stack)
                stack=stack(:,:,1:zz);
            end;
            break;
        end;
        zz=zz+1;
        
    end;

    return;
end;






    