function tom_av2_build_file_struct(directory,flag,number)
%TOM_AV2_BUILD_FILE_STRUCT creates ...
%
%   tom_av2_build_file_struct(directory,flag,number)
%
%PARAMETERS
%
%  INPUT
%   directory           ...
%   flag                ...
%   number              ...
%  
%  OUTPUT
%
%EXAMPLE
%   tom_av2_build_file_struct(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   ...
%
%   created by .. mm/dd/yy
%   updated by ..
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
    flag='all';
 end;

if nargin < 1
    directory='./';
end;
    
   
if nargin < 3
    number=5;
 end;

    
 
%rock the struct


        
        
 if  (strcmp(flag,'align') | strcmp(flag,'all'))  
        
    mkdir align/low/;
        mkdir align/low/whole_particle/;
        mkdir align/low/core/;
     mkdir align;
        mkdir align/high/;
        mkdir align/high/whole_particle/;
        mkdir align/high/core/;
 
 end;
        
if  (strcmp(flag,'pick') | strcmp(flag,'all'))    
        
    mkdir pick;
    mkdir pick/high/;
    mkdir pick/low/;
end;       

    
if  (strcmp(flag,'org') | strcmp(flag,'all'))    

    mkdir org;
end;


if  (strcmp(flag,'sort') | strcmp(flag,'all'))    

    mkdir sort;
    mkdir sort/high/;
    mkdir sort/low/; 
end;
    

if  (strcmp(flag,'rec') | strcmp(flag,'all'))    
            
        if (exist([directory '/logs'],'dir')==0);
            mkdir([directory '/logs']);
        end;
    
    
        for i=1:number
            if exist([directory '/step' num2str(i)],'dir')==0;
                 mkdir([directory '/step' num2str(i)]);
            end;
            if exist([directory '/step' num2str(i) '/proj'],'dir')==0;
                mkdir([directory '/step' num2str(i) '/proj' ]);
            end;
            
            if exist([directory '/step' num2str(i) '/avg'],'dir')==0; 
                mkdir([directory '/step' num2str(i) '/avg' ]);
            end;
            
            if exist([directory '/step' num2str(i) '/model'],'dir')==0; 
                mkdir([directory '/step' num2str(i) '/model' ]);
            end;
            
             if exist([directory '/step' num2str(i) '/align'],'dir')==0; 
                mkdir([directory '/step' num2str(i) '/align' ]);
            end;
            
        end;
end;
        



if  (strcmp(flag,'rec-logs') )
    if (exist([directory '/logs'],'dir')==0);
            mkdir([directory '/logs']);
    end;
end;    
    




