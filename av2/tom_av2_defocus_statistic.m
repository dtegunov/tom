function st=tom_av2_defocus_statistic(path,filename,extension)


disp('scanning files ...');

[dircell, flagcell_out] = get_dircontents([path '/' filename '*' extension],{},{},path);

%build info structure 
num_of_clusters=0;
for i=1:size(dircell,2)
    h=tom_reademheader(dircell{i});
    [tmp] = parse_header(h.Header,dircell{i});
    info(i)=tmp;
    if (info(i).cluster > num_of_clusters)
        num_of_clusters=info(i).cluster;
    end;

end;


for i=1:size(info,2)
    def(i)=info(i).fitted_def;
    if (def(i)==1)
        disp(info(i).filename);
    end;
end;

[a b c d]=tom_dev(def);
st.mean=a;
st.min=b;
st.max=c;
st.std=d;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  parse header Informatin                            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [info] = parse_header(header,filename)



[xx,pos] = strtok(char(header.Comment), ';');
[pos,remain] = strtok(pos,';');
[xx,nominal] = strtok(remain,';');
[nominal,cluster] = strtok(nominal,'!');
nominal = nominal(2:end);
cluster = cluster(2:end);
[pos_x,remain] = strtok(pos,' ');
[pos_y,pos_z] = strtok(remain,' ');

 fitted_def=header.FocusIncrement;

 info.pos(1)=str2num(pos_x);
 info.pos(2)=str2num(pos_y);
 info.pos(3)=str2num(pos_z);
 info.mesured_def=str2num(nominal);
 info.intended_def=header.Defocus./10000;
 info.fitted_def=fitted_def./10000;
 info.cluster=str2num(cluster);
 info.filename=filename;
 if (header.Marker_X==999)
    info.isinterp=1;
 else
    info.isinterp=0; 
 end;
 
 if (fitted_def~=1)
     info.isfitted=1;
 else
     info.isfitted=0;
 end;

 
if isempty(findstr(nominal,';'))==0
    info.mesured_def=-5;
    info.fitted_def=1;
    info.isfitted=0;
end;

if isempty(info.mesured_def)
    info.isfitted=0;
    info.mesured_def=-5;
    info.fitted_def=1;
end;


if (info.isinterp==1)
    info.isfitted=0;
    %info.mesured_def=1;
    info.fitted_def=1;
end;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  get all the em-files in a directory                                %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dircell, flagcell_out] = get_dircontents(directory, dircell, flagcell,path)

%Get list of EM-Files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%h = waitbar(0,'Getting directory list...');

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
            if tom_isemfile([path '/' dirlist(i).name]) == 1
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
    %waitbar(i./files,h,[num2str(i), ' of ', num2str(files), ' files scanned.']);
end

%close(h);

if size(dircell,2) == 0
    errordlg('No EM files could be found in this directory.');
    return;
end

%set(findobj('Tag','fileslider'),'Max',size(dircell,2),'SliderStep',[1./size(dircell,2) 1./size(dircell,2)]);
