function vect=tom_exract_num_from_file(filename)
% vect=tom_exract_num_from_file(filename) gives vector with indices back
%
%    vect=tom_exract_num_from_file(filename)
%
%   tom_exract_num_from_file gives vector with indices back
%  
%  
%
%PARAMETERS
%
%  INPUT
%   filename           input filename
%   
%
%
%  OUTPUT
%    vect           indices vect     
%    
%
%
%EXAMPLE
%     



in_f=importdata(filename);

vect=zeros(1,length(in_f.textdata));
for i=1:length(in_f.textdata)
    % create high listlines{zz}=b;
    [a b c]=fileparts(in_f.textdata{i});
    
    if (findstr(b,'_')==1)
        [a d]=strtok(b,'_');
        vect(i)=str2double(strrep(d,'_',''));
    else
        for ii=1:length(b)
            if (isnan(str2double(b(ii)))==0 )
                vect(i)=str2double(b(ii:end));
                break;
            end
        end;
    end;
        
end;