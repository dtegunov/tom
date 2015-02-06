function tom_av2_ctf_interpol(path,filename,extension)
%tom_av2_ctf_interpol interpolates the focus of given dataset
%the projections which should be interpolated have to have 1
%header.FocusIncrement. Interpolated projections are taged with 999
%in header.MarkerX. 
%
%
%tom_av2_ctf_interpol(path,filename,extension)
%
%PARAMETERS
%
%  INPUT
%   path               path of the projections
%   filename           filename of the projections              
%   extension          extension of the projections  
%
%
%
%EXAMPLE
%
%tom_av2_ctf_interpol('./','26S','.em');
%
%
%REFERENCES
%
%SEE ALSO
%   tom_ctftool
%
%   created by fb 31/01/07
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

disp('calculating surfaces ...');

%calculate surface
sum_fit_def=0;
fit_count=1;
for ii=1:num_of_clusters
    zz=1;
    zz2=1;
    for i=1:size(info,2)
        if (info(i).cluster==ii)
            mesh_x(zz2)=info(i).pos(1);
            mesh_y(zz2)=info(i).pos(2);
            zz2=zz2+1;
            if (info(i).isfitted==1)
                data_x(zz)=info(i).pos(1);
                data_y(zz)=info(i).pos(2);
                data_z(zz)=(info(i).fitted_def-info(i).intended_def ) + (info(i).pos(3) + info(i).mesured_def);
                sum_fit_def(fit_count)=info(i).fitted_def;
                fit_count=fit_count+1;
                zz=zz+1;
            end;
        end;
    end;
    if ~isempty(data_x)
        
        [X,Y]=meshgrid(sort(mesh_x),sort(mesh_y));
        warning('off'); zi{ii} = griddata(data_x,data_y,data_z,X,Y); warning('on');
        Xi{ii}=round(X.*10000);
        Yi{ii}=round(Y.*10000);
    end
    zz=1;
    zz2=1;
    mesh_x = [];
    mesh_y = [];
    data_x = [];
    data_y = [];
    data_z = [];
%     clear('mesh_x');
%     clear('mesh_y');
%     clear('data_x');
%     clear('data_y');
%     clear('data_z');
    disp(num2str(ii));
end;

disp(['Fitting Statistics: ' num2str(size(sum_fit_def,2) ) ' of ' num2str(size(info,2)) ' fitted']);
mean_fitted=tom_dev(sum_fit_def);


disp('interpolating ');

zz=1;

% calclate iterpolated defocus
for i=1:size(info,2)
    
    clu=info(i).cluster;
    if (info(i).isfitted==0)
        if isempty(clu)
            idx=[];
        else
            try;idx=find( (Xi{clu}==(info(i).pos(1)*10000))  & (Yi{clu}==info(i).pos(2).*10000) );catch idx = [];end
        end;
        
        if (isempty(idx)==0)
            tmp=zi{clu};
            interp=tmp(idx(1));
            int_fit_def=interp-info(i).mesured_def-info(i).pos(3)+info(i).intended_def;
        else
            int_fit_def=mean_fitted;
        end;
        
        if int_fit_def < (mean_fitted.*1.3) | int_fit_def > (mean_fitted./1.3) |  isnan(int_fit_def)
            int_fit_def=mean_fitted;
        end;
        
        %updating Header
        h=tom_reademheader(info(i).filename);
        h.Header.FocusIncrement=int_fit_def.*10000;
        h.Header.Marker_X=999;
        tom_writeemheader(info(i).filename,h.Header);
        
        tmp_interp_def(zz)=int_fit_def;
        zz=zz+1;
        
    end;
    
end;

disp('end');



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

set(findobj('Tag','fileslider'),'Max',size(dircell,2),'SliderStep',[1./size(dircell,2) 1./size(dircell,2)]);
