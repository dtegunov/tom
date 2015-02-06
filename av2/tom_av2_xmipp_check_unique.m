function [isunique lines]=tom_av2_xmipp_check_unique(filename)
%tom_av2_xmipp_check_unique(filename) checks if entris are unique
%
%   tom_av2_xmipp_check_unique(filename)
%
%   tom_av2_xmipp_check_unique checks if partnr in a sel doc or htl file
%   are unique ...use if xmipp is crashing with that error
%  
%
%PARAMETERS
%
%  INPUT
%   input_sel           sel,doc or htl fiel
%   
%
%
%  OUTPUT
%    isunique          0 for nonunique 1 for unique
%    lines             the grept files   
%
%
%EXAMPLE
%     
%  
% NOTE:
%  
% The files have to have the right extension ... reading is triggert!!
% 
%
%
%


isunique=1;
lines='';

[a b c]=fileparts(filename);

if ((strcmp(c,'.htl') || strcmp(c,'.sel') || strcmp(c,'.doc'))==0 )
    error('Wrong data format only htl,sel and doc supported');
end;


disp(' ');
disp('...POWER CHECKER ROCKING ');

disp(' ');

if (strcmp(c,'.htl'))
    
    disp(['reading ' filename]);
    in_htl=importdata(filename);
    disp(' ');
    disp(['coverting 2 index ' filename]);
    disp(' ');
    
    %converting htl 2 index
    names_h_tmp=zeros(length(in_htl),1);
    names_l_tmp=zeros(length(in_htl),1);
    
    for i=1:length(in_htl)
        % create high listlines{zz}=b;
        
        [start rest]=strtok(in_htl{i});
        
        names_h{i}=deblank(start);
        [a b c]=fileparts(names_h{i});
        [a d]=strtok(b,'_');
        names_h_tmp(i)=str2double(strrep(d,'_',''));
        
        % create low list
        tmp_l=deblank(strrep(rest,',',''));
        names_l{i}=tmp_l(2:end);
        [a b c]=fileparts(names_l{i});
        [a d]=strtok(b,'_');
        names_l_tmp(i)=str2double(strrep(d,'_',''));
     end;
     
     
     if (length(unique(names_l_tmp) ) ==i)
         disp('High column part numbers in htl are unique !');
     else
         disp('ERROR: low column is not unique!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
         disp(['length uniuqe: ' num2str(length(unique(names_l_tmp))) ' length totoal: ' num2str(length(names_l_tmp)) ]);
         disp(' ');
         isunique=0;
         [tmp_sort]=sort(names_l_tmp);
         idx=find(diff(tmp_sort)==0);
         zz=1;
         for ii=1:length(idx)
            call=['grep -n ' num2str(tmp_sort(idx(ii))) ' ' filename];
            [a b]=unix(call);
            lines{zz}=b;
            disp(lines{zz});
            zz=zz+1;
           disp(' ');
         end;
     end;
     
     if (length(unique(names_h_tmp) ) ==i)
         disp('Low column part numbers in htl are unique !');
     else
         disp('ERROR: high column is not unique!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
         disp(['length uniuqe: ' length(unique(names_h_tmp)) ' length totoal: ' length(names_h_tmp)]);
         isunique=0;
         [tmp_sort]=sort(names_h_tmp);
         idx=find(diff(tmp_sort)==0);
         zz=1;
         for ii=1:length(idx)
            call=['grep -n ' num2str(tmp_sort(idx(ii))) ' ' filename];
            [a b]=unix(call);
            lines{zz}=b;
            disp(lines{zz});
            zz=zz+1;
           disp(' ');
         end;  
    end;
     
end;



if (strcmp(c,'.sel'))
    disp(['reading ' filename]);
    in_sel=importdata(filename);
    disp(['coverting 2 index ' filename]);
    
    %converting input sel 2 index
    in_sel_tmp=zeros(length(in_sel.textdata),1);
    for ii=1:length(in_sel.textdata)
        [a b c]=fileparts(in_sel.textdata{ii});
        [a d]=strtok(b,'_');
        in_sel_tmp(ii)=str2double(strrep(d,'_',''));
    end;
    
    disp(' ');
    if (length(unique(in_sel_tmp) ) ==ii)
        disp('part numbers in sel are unique !');
    else
         disp('ERROR: sel is not unique!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
         disp(['length uniuqe: ' num2str(length(unique(in_sel_tmp))) ' length totoal: ' num2str(length(in_sel_tmp)) ]);
         disp(' ');
         isunique=0;
         [tmp_sort]=sort(in_sel_tmp);
         idx=find(diff(tmp_sort)==0);
         zz=1;
         for ii=1:length(idx)
            call=['grep -n ' num2str(tmp_sort(idx(ii))) ' ' filename];
            [a b]=unix(call);
            lines{zz}=b;
            disp(lines{zz});
            zz=zz+1;
           disp(' ');
         end;
    end;
    
end;



if (strcmp(c,'.doc'))
    disp('Not implemented !!!!!!!!!!');
    
end;















