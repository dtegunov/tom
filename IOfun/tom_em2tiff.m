function tom_em2tiff(path_in,path_out,number_of_files)
%TOM_EM2TIFF creates ...
%
%   tom_em2tiff(path_in,path_out,number_of_files)
%
%   Converts micrographs stored in EM format to TIFF
%
%PARAMETERS
%
%  INPUT
%   optional parameters are :
%
%   path_in             ...
%   path_out            ...
%   number_of_files     ...
%  
%  OUTPUT
%
%EXAMPLE
%   tom_em2tiff();
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by ... (author date)
%   updated by TH2 (21.12.09)
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

% convert EM format files to TIF format

if nargin == 3
    for i=1:number_of_files
        name_in=[path_in num2str(i) '.em' ];
        name_out=[path_out num2str(i) '.tif'];
        i=tom_emread(name_in);
        i=tom_norm(i.Value,1);
        imwrite(i',name_out,'tiff');
    end;
elseif nargin == 0
    error('Type "help tom_em2tiff" for how to use this function!');
else    
    [files path]= uigetfile({'*.em'},'MultiSelect', 'on','Select one or more em files in one folder.');
    
    
    if (~iscell(files) && ~isstr(files))
        error('You must select a file!');
    end;
    
    if ~iscell(files)
        files = {files};
    end;   
    
    destination = uigetdir('Select destination directory');
    
    seperator = '/';
    
    if ispc
        seperator = '\';
    end;
    
    for filesIterator = 1:length(files)
       file = files{filesIterator};
       
       outName = file(1:end-3);
       
       image=tom_emread([path file]);
       image=tom_norm(image.Value,1);
       
       imwrite(image',[destination seperator outName '.tiff'],'tiff');
      
    end;
end;
