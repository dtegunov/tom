function tom_adjust_positions(logfile,directory)

%read log file
fid = fopen(logfile);

linedata = 1;
zeile = 1;
while 1

    %read a line from the file
    tline = fgetl(fid);
    if tline == -1
        break;
    end
    
    %search for 11 colons
    idx = strfind(tline,';');
    if isempty(idx) || length(idx) < 11
        continue;
    end
    
    %extract data
    lauf = 1;
    while true
       [str, tline] = strtok(tline, ';');
       if isempty(str),  break;  end
       if ~isempty(str2num(deblank(str)))
           linedata(zeile,lauf) =  str2num(deblank(str)); 
           lauf=lauf+1;
       end
    end
    zeile = zeile+1;
end

fclose(fid);

% pos_x = [];
% pos_y = [];
% pos_z = [];
% 
% pos_x = linedata(:,6);
% pos_y = linedata(:,7);
% pos_z = linedata(:,8);
% 
pos_new = linedata(:,6:8);
point = linedata(:,2);
clear('linedata');

%correct all comments
[dircell] = get_dircontents(directory, {}, {});
h = waitbar(0,'Getting directory list...');
for i=1:length(dircell)
    file = dircell{i};
    header = tom_reademheader(file);
    comment = header.Header.Comment;
    rem = deblank(char(comment));
    
    for i=1:3;[s{i},rem]=strtok(rem,';');end
    pos_old = str2num(s{2});
    idx = tom_nearestpoint(pos_old,pos_new);
    
    s{2} = [num2str(pos_new(idx,1)) ' ' num2str(pos_new(idx,2)) ' ' num2str(pos_new(idx,3))];
    
    header.Header.Comment = [s{1} ';' s{2} ';' s{3} deblank(rem) num2str(point(idx))];
    tom_writeemheader(file,header.Header);
    waitbar(i./length(dircell),h,[num2str(i), ' of ', num2str(length(dircell)), ' files corrected.']);
end
close(h);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  get all the em-files in a directory                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dircell, flagcell_out] = get_dircontents(directory, dircell, flagcell)

%Get list of EM-Files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = waitbar(0,'Getting directory list...');

dirlist = dir(directory);
files = size(dirlist,1);
j = size(dircell,2);
if j == 0
    j = 1;
end

%sort out all the em-files
for i = 1:files
    if dirlist(i).isdir == 0
        if isempty(strmatch(dirlist(i).name,dircell))
            if tom_isemfile([directory '/' dirlist(i).name]) == 1
                dircell{j} = dirlist(i).name;
                try
                    flagcell_out{j} = flagcell{j};
                catch
                    flagcell_out{j} = 1;
                end
                j = size(dircell,2) + 1;
            end
        end
    end
    waitbar(i./files,h,[num2str(i), ' of ', num2str(files), ' files scanned.']);
end

close(h);

if size(dircell,2) == 0
    errordlg('No EM files could be found in this directory.');
    return;
end

set(findobj('Tag','fileslider'),'Max',size(dircell,2),'SliderStep',[1./size(dircell,2) 1./size(dircell,2)]);
