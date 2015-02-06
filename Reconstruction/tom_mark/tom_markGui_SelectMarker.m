function varargout = tom_markGui_SelectMarker(varargin)
% No help (yet?)


gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @figure_OpeningFcn, ...
                   'gui_OutputFcn',  @figure_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% --- Executes just before tom_markGui_SelectMarker is made visible.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function figure_OpeningFcn(hObject, eventdata, handles, varargin)

markGui_hfigure = findall(0, 'Type', 'figure', 'Tag', 'figure_markGui');
if (length(markGui_hfigure) ~= 1)
    %errordlg('Call tom_markGui_Alignment3dItr from tom_markGui', 'markGui_Alignment3dItr', 'modal');
    warning('Call tom_markGui_SelectMarker from tom_markGui');
    tom_markGui;
    return;
end;
markGui_data = guidata(markGui_hfigure);


if (length(markGui_data.data.filenames) < 2)
    error('No images loaded in tom_markGui');
end;
if (isempty(markGui_data.data.markersets) || isempty(markGui_data.data.markersets(markGui_data.data.sel_markerset).markerset))
    error('No markerset loaded in tom_markGui');
end;

[tiltangles, imsize] = markGui_data.fcn_local.getTiltangles(markGui_data.data.emheader);


if (~isfield(markGui_data.data,'markGui_SelectMarker') || ~isstruct(markGui_data.data.markGui_SelectMarker) || numel(markGui_data.data.markGui_SelectMarker)~=1)
    handles.data.cfg = struct();
else
    handles.data.cfg = markGui_data.data.markGui_SelectMarker;
end;
    

if (~isfield(handles.data.cfg, 'seltype_replace') || ~islogical(handles.data.cfg.seltype_replace) || numel(handles.data.cfg.seltype_replace)~=1)
    handles.data.cfg.seltype_replace = false;
end;


if (~isfield(handles.data.cfg, 'nlength_min') || ~islogical(handles.data.cfg.nlength_min) || numel(handles.data.cfg.nlength_min)~=1)
    handles.data.cfg.nlength_min = true;
end;
if (~isfield(handles.data.cfg, 'nlength_minval') || ~isnumeric(handles.data.cfg.nlength_minval) || numel(handles.data.cfg.nlength_minval)~=1)
    handles.data.cfg.nlength_minval = 1;
end;
if (~isfield(handles.data.cfg, 'nlength_max') || ~islogical(handles.data.cfg.nlength_max) || numel(handles.data.cfg.nlength_max)~=1)
    handles.data.cfg.nlength_max = false;
end;
if (~isfield(handles.data.cfg, 'nlength_maxval') || ~isnumeric(handles.data.cfg.nlength_maxval) || numel(handles.data.cfg.nlength_maxval)~=1)
    handles.data.cfg.nlength_maxval = length(markGui_data.data.filenames);
end;
if (handles.data.cfg.nlength_maxval > length(markGui_data.data.filenames))
    handles.data.cfg.nlength_maxval = length(markGui_data.data.filenames);
end;
if (handles.data.cfg.nlength_minval > handles.data.cfg.nlength_maxval)
    handles.data.cfg.nlength_minval = handles.data.cfg.nlength_maxval;
end;
if (~isfield(handles.data.cfg, 'align_markerset_idx') || ~isnumeric(handles.data.cfg.align_markerset_idx) || numel(handles.data.cfg.align_markerset_idx)~=1 || handles.data.cfg.align_markerset_idx~=round(handles.data.cfg.align_markerset_idx) || ~(handles.data.cfg.align_markerset_idx>0))
    handles.data.cfg.align_markerset_idx = 1;
end;
if (~isfield(handles.data.cfg, 'align_markerset_name') || ~ischar(handles.data.cfg.align_markerset_name))
    handles.data.cfg.align_markerset_name = '';
end;
if (~isfield(handles.data.cfg, 'align_markerset_alignment') || ~ischar(handles.data.cfg.align_markerset_alignment))
    handles.data.cfg.align_markerset_alignment = '';
end;
if (~isfield(handles.data.cfg, 'align_dist_mean') || ~islogical(handles.data.cfg.align_dist_mean) || numel(handles.data.cfg.align_dist_mean)~=1)
    handles.data.cfg.align_dist_mean = false;
end;
if (~isfield(handles.data.cfg, 'align_dist_meanval') || ~isnumeric(handles.data.cfg.align_dist_meanval) || numel(handles.data.cfg.align_dist_meanval)~=1)
    handles.data.cfg.align_dist_meanval = 50;
end;
if (~isfield(handles.data.cfg, 'align_dist_max') || ~islogical(handles.data.cfg.align_dist_max) || numel(handles.data.cfg.align_dist_max)~=1)
    handles.data.cfg.align_dist_max = false;
end;
if (~isfield(handles.data.cfg, 'align_dist_maxval') || ~isnumeric(handles.data.cfg.align_dist_maxval) || numel(handles.data.cfg.align_dist_maxval)~=1)
    handles.data.cfg.align_dist_maxval = 50;
end;
if (~isfield(handles.data.cfg, 'align_dist3d') || ~islogical(handles.data.cfg.align_dist3d) || numel(handles.data.cfg.align_dist3d)~=1)
    handles.data.cfg.align_dist3d = false;
end;
if (~isfield(handles.data.cfg, 'align_dist3dval') || ~isnumeric(handles.data.cfg.align_dist3dval) || numel(handles.data.cfg.align_dist3dval)~=1)
    handles.data.cfg.align_dist3dval = mean(imsize(:))/ 2;
end;
if (~isfield(handles.data.cfg, 'align_dist3d_regplane') || ~islogical(handles.data.cfg.align_dist3d_regplane) || numel(handles.data.cfg.align_dist3d_regplane)~=1)
    handles.data.cfg.align_dist3d_regplane = false;
end;
if (~isfield(handles.data.cfg, 'align_dist3d_regplaneval') || ~isnumeric(handles.data.cfg.align_dist3d_regplaneval) || numel(handles.data.cfg.align_dist3d_regplaneval)~=1)
    handles.data.cfg.align_dist3d_regplaneval = mean(imsize(:))/ 4;
end;




set(handles.radiobutton_seltype_add, 'Value', ~handles.data.cfg.seltype_replace);
set(handles.radiobutton_seltype_replace, 'Value', handles.data.cfg.seltype_replace);

set(handles.checkbox_nlength_min, 'Value', handles.data.cfg.nlength_min);
set(handles.edit_nlength_min, 'String', num2str(handles.data.cfg.nlength_minval));
set(handles.checkbox_nlength_max, 'Value', handles.data.cfg.nlength_max);
set(handles.edit_nlength_max, 'String', num2str(handles.data.cfg.nlength_maxval));

set(handles.checkbox_align_dist_max, 'Value', handles.data.cfg.align_dist_max);
set(handles.edit_align_dist_max, 'String', num2str(handles.data.cfg.align_dist_maxval));
set(handles.checkbox_align_dist_mean, 'Value', handles.data.cfg.align_dist_mean);
set(handles.edit_align_dist_mean, 'String', num2str(handles.data.cfg.align_dist_meanval));
set(handles.checkbox_align_dist3d, 'Value', handles.data.cfg.align_dist3d);
set(handles.edit_align_dist3d, 'String', num2str(handles.data.cfg.align_dist3dval));
set(handles.checkbox_align_dist3d_regplane, 'Value', handles.data.cfg.align_dist3d_regplane);
set(handles.edit_align_dist3d_regplane, 'String', num2str(handles.data.cfg.align_dist3d_regplaneval));

handles.data.markGui_data = markGui_data;


[handles.data.align.markerset_names, handles.data.align.markerset_alignfields, handles.data.align.markerset_alignfields_text, handles.data.align.markerset_alignfields_valid] = getAlignments(handles);


if (length(handles.data.align.markerset_names) >= handles.data.cfg.align_markerset_idx && strcmp(handles.data.align.markerset_names{handles.data.cfg.align_markerset_idx}, handles.data.cfg.align_markerset_name))
    sel_align_markerset = handles.data.cfg.align_markerset_idx;
else
    sel_align_markerset = handles.data.markGui_data.data.sel_markerset;
    handles.data.cfg.align_alignment = '';
end;

sel = find(handles.data.align.markerset_alignfields_valid == sel_align_markerset);
if (isempty(sel))
    sel = 1;
    sel_align_markerset = [];
end;


set(handles.listbox_align_markersets, 'String', handles.data.align.markerset_names(handles.data.align.markerset_alignfields_valid), 'Value', sel);

Callback_listbox_align_markersets(handles.listbox_align_markersets, [], handles);

if (length(sel_align_markerset)==1)
    sel_alignment = find(strcmp(handles.data.align.markerset_alignfields{sel_align_markerset}, handles.data.cfg.align_markerset_alignment));
    if (length(sel_alignment) == 1)
        set(handles.listbox_align_alignments, 'Value', sel_alignment);
    end;
    enable_align = 'on';
else
    enable_align = 'off';
end;
set(handles.checkbox_align_dist_max, 'Enable', enable_align);
set(handles.checkbox_align_dist_mean, 'Enable', enable_align);
set(handles.edit_align_dist_max, 'Enable', enable_align);
set(handles.edit_align_dist_mean, 'Enable', enable_align);
set(handles.pushbutton_align_dist_max, 'Enable', enable_align);
set(handles.pushbutton_align_dist_mean, 'Enable', enable_align);

set(handles.checkbox_align_dist3d, 'Enable', enable_align);
set(handles.edit_align_dist3d, 'Enable', enable_align);
set(handles.pushbutton_align_dist3d, 'Enable', enable_align);
set(handles.checkbox_align_dist3d_regplane, 'Enable', enable_align);
set(handles.edit_align_dist3d_regplane, 'Enable', enable_align);
set(handles.pushbutton_align_dist3d_regplane, 'Enable', enable_align);
set(handles.pushbutton_align_dist3d_regplane_view, 'Enable', enable_align);


guidata(hObject, handles);



uiwait(handles.figure_markGui_SelectMarker);


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% --- Outputs from this function are returned to the command line.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = figure_OutputFcn(hObject, eventdata, handles) 


varargout = cell(1,2);

if (isempty(handles))
    return;
end;

if (isfield(handles,'output'))
    varargout{1} = handles.output.cfg;
    varargout{2} = handles.output.sel_markers;
end;
if (ishandle(handles.figure_markGui_SelectMarker))
    delete(handles.figure_markGui_SelectMarker); 
end;


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Callback_listbox_align_markersets(hObject, eventdata, handles)

sel = get(handles.listbox_align_markersets, 'Value');

if (isempty(sel) || isempty(get(handles.listbox_align_markersets,'String')) || length(handles.data.align.markerset_alignfields_valid)<sel)
    set(handles.listbox_align_alignments, 'String', {}, 'Value', 1);
else
    set(handles.listbox_align_alignments, 'String', handles.data.align.markerset_alignfields_text{handles.data.align.markerset_alignfields_valid(sel)}, 'Value', 1);
end;
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Callback_pushbutton_view(hObject, eventdata, handles, type)

ftag = ['tom_markGui_SelectMarker_view'];
h = findobj('Type', 'figure', 'Tag', ftag);
if (isempty(h))
    h = figure();
    set(h, 'Tag', ftag);
end;
figure(h);
cla;


switch (type)
    case {'dist3d_regplane_view'}
        [align_markerset_idx, align_markerset_name, align_markerset_alignment, aligncnf] = getSelectedAlignment(handles);
        if (isempty(aligncnf))
            delete(h);
            return;
        end;
        
        X = aligncnf.X(1:3,:)';
        X = X(all(isfinite(X),2),:);
        % http://www.mathworks.com/products/statistics/demos.html?file=/products/demos/shipping/stats/orthoregdemo.html
        [coeff,score,roots] = princomp(X);
        basis = coeff(:,1:2);
        normal = coeff(:,3);
        pctExplained = roots' ./ sum(roots);
 
        [n,p] = size(X);
        meanX = mean(X,1);
        Xfit = repmat(meanX,n,1) + score(:,1:2)*coeff(:,1:2)';
        residuals = X - Xfit;
        
        error = abs((X - repmat(meanX,n,1))*normal);
        sse = sum(error.^2);
        
        [xgrid,ygrid] = meshgrid(linspace(min(X(:,1)),max(X(:,1)),5), linspace(min(X(:,2)),max(X(:,2)),5));
        zgrid = (1/normal(3)) .* (meanX*normal - (xgrid.*normal(1) + ygrid.*normal(2)));
        h = mesh(xgrid,ygrid,zgrid,'EdgeColor',[0 0 0],'FaceAlpha',0);

        hold on
        above = (X-repmat(meanX,n,1))*normal > 0;
        below = ~above;
        nabove = sum(above);
        X1 = [X(above,1) Xfit(above,1) nan*ones(nabove,1)];
        X2 = [X(above,2) Xfit(above,2) nan*ones(nabove,1)];
        X3 = [X(above,3) Xfit(above,3) nan*ones(nabove,1)];
        plot3(X1',X2',X3','-', X(above,1),X(above,2),X(above,3),'o', 'Color',[0 .7 0
        ]);
        nbelow = sum(below);
        X1 = [X(below,1) Xfit(below,1) nan*ones(nbelow,1)];
        X2 = [X(below,2) Xfit(below,2) nan*ones(nbelow,1)];
        X3 = [X(below,3) Xfit(below,3) nan*ones(nbelow,1)];
        plot3(X1',X2',X3','-', X(below,1),X(below,2),X(below,3),'o', 'Color',[1 0 0]);

        hold off
        maxlim = max(abs(X(:)))*1.1;
        axis([-maxlim maxlim -maxlim maxlim -maxlim maxlim]);
        axis equal
        view(-23.5,5);
        
    case {'distance_mean', 'distance_max', 'dist3d', 'dist3d_regplane'}
        [align_markerset_idx, align_markerset_name, align_markerset_alignment, aligncnf] = getSelectedAlignment(handles);
        if (isempty(aligncnf))
            delete(h);
            return;
        end;
        
        switch (type)
            case {'distance_mean'}
                vals = aligncnf.distance_mean;
                n = str2double(get(handles.edit_align_dist_mean,'String'));
            case {'distance_max'}
                vals = aligncnf.distance_max;
                n = str2double(get(handles.edit_align_dist_max,'String'));
            case {'dist3d'}
                X = aligncnf.X(1:3,:);
                X = X(:,all(isfinite(X),1));
                vals = sqrt(sum((X(1:3,:) - repmat(mean(X(1:3,:),2), [1, size(X,2)])) .^2, 1));
                n = str2double(get(handles.edit_align_dist3d,'String'));
            case {'dist3d_regplane'}
                X = aligncnf.X(1:3,:);
                X = X(:,all(isfinite(X),1));
                
                p = princomp(X(1:3,:)');
                meanX = mean(X(1:3,:),2);
                p = p(:,3);
                p = [p(1:3); - p'*meanX];
                
                vals = abs(p(1:3)' * (X(1:3,:) - repmat(meanX, [1, size(X,2)])));
                n = str2double(get(handles.edit_align_dist3d_regplane,'String'));
        end;     
        vals = sort(vals);
        
        tom_dev(vals(isfinite(vals)));
        %vals = sort(vals, 'ascend');
        hold on;
        plot(1:length(vals), vals, 'b-');
        if (numel(n)==1 && isfinite(n))
            plot([1,length(vals)], [n,n],'r-');
        end;
        hold off;
        
        
    case 'nlength'
        markerset = handles.data.markGui_data.data.markersets(handles.data.markGui_data.data.sel_markerset).markerset;
        nlength = squeeze(sum(all(isfinite(markerset(1:2,:,:)), 1), 2))';
        tom_dev(nlength);
        
        edges = 0:1:size(markerset, 2);
        nlengthc = histc(nlength, edges);
        bar(edges, nlengthc, 'histc');
    otherwise
        warning(['Unrecognized type ' type]);
        delete(h);
        return;
end;

 








end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ss, nother, notherverbose, useidx] = getAlignments(handles)


nother = cell(1, length(handles.data.markGui_data.data.markersets));
ss = cell(1, length(nother));
useidx = [];
for (i=1:length(nother))
    notheri = {};
    if (isfield(handles.data.markGui_data.data.markersets(i), 'alignment3d') && ...
        isstruct(handles.data.markGui_data.data.markersets(i).alignment3d) && ...
        numel(handles.data.markGui_data.data.markersets(i).alignment3d)==1 && ...
        isfield(handles.data.markGui_data.data.markersets(i).alignment3d, 'done') && ...
        handles.data.markGui_data.data.markersets(i).alignment3d.done)
        otheralignments = true;
        notheri = { notheri{:}, 'alignment3d'};
    end;
    if (isfield(handles.data.markGui_data.data.markersets(i), 'alignment3ditr') && ...
        isstruct(handles.data.markGui_data.data.markersets(i).alignment3ditr) && ...
        numel(handles.data.markGui_data.data.markersets(i).alignment3ditr)==1 && ...
        isfield(handles.data.markGui_data.data.markersets(i).alignment3ditr, 'done') && ...
        handles.data.markGui_data.data.markersets(i).alignment3ditr.done)
        notheri = { notheri{:}, 'alignment3ditr'};
    end;
    nother{i} = notheri;
    if (~isempty(notheri))
        useidx(end+1) = i;
    end;
    ss{i} = [num2str(i) ': ' handles.data.markGui_data.data.markersets(i).name ' (' num2str(size(handles.data.markGui_data.data.markersets(i).markerset,3)) ' marker)'];
end;

notherverbose = cell(size(nother));
for (i=1:length(nother))
    notheri = nother{i};
    notherverbosei = {};
    for (j=1:length(notheri))
        switch (notheri{j})
            case 'alignment3d'
                notherverbosei{j} = 'Rigid Body Alignmentd3d';
            case 'alignment3ditr'
                notherverbosei{j} = 'Rigid Body Alignmentd3d (iterative)';
            otherwise
                notherverbosei{j} = notheri{j};
        end;        
    end; 
    notherverbose{i} = notherverbosei;
end;



end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Callback_pushbutton_ok(hObject, eventdata, handles)

cfg = handles.data.cfg;

cfg.seltype_replace = logical(get(handles.radiobutton_seltype_replace, 'Value'));

cfg.nlength_min = logical(get(handles.checkbox_nlength_min, 'Value'));
cfg.nlength_max = logical(get(handles.checkbox_nlength_max, 'Value'));
if (cfg.nlength_min)
    n = str2double(get(handles.edit_nlength_min, 'String'));
    if (isfinite(n))
        cfg.nlength_minval = n;
    else
        cfg.nlength_min = false;
    end;
end;
if (cfg.nlength_max)
    n = str2double(get(handles.edit_nlength_max, 'String'));
    if (isfinite(n))
        cfg.nlength_maxval = n;
    else
        cfg.nlength_max = false;
    end;
end;
cfg.align_dist_max = logical(get(handles.checkbox_align_dist_max,'Value'));
cfg.align_dist_mean = logical(get(handles.checkbox_align_dist_mean,'Value'));
cfg.align_dist3d = logical(get(handles.checkbox_align_dist3d,'Value'));
cfg.align_dist3d_regplane = logical(get(handles.checkbox_align_dist3d_regplane,'Value'));
if (cfg.align_dist_max)
    n = str2double(get(handles.edit_align_dist_max, 'String'));
    if (isfinite(n))
        cfg.align_dist_maxval = n;
    else
        cfg.align_dist_max = false;
    end;
end;
if (cfg.align_dist_mean)
    n = str2double(get(handles.edit_align_dist_mean, 'String'));
    if (isfinite(n))
        cfg.align_dist_meanval = n;
    else
        cfg.align_dist_mean = false;
    end;
end;
if (cfg.align_dist3d)
    n = str2double(get(handles.edit_align_dist3d, 'String'));
    if (isfinite(n))
        cfg.align_dist3dval = n;
    else
        cfg.align_dist3d = false;
    end;
end;
if (cfg.align_dist3d_regplane)
    n = str2double(get(handles.edit_align_dist3d_regplane, 'String'));
    if (isfinite(n))
        cfg.align_dist3d_regplaneval = n;
    else
        cfg.align_dist3d_regplane = false;
    end;
end;

[cfg.align_markerset_idx, cfg.align_markerset_name, cfg.align_markerset_alignment, aligncnf] = getSelectedAlignment(handles);

markerset = handles.data.markGui_data.data.markersets(handles.data.markGui_data.data.sel_markerset).markerset;

sel_markers = true(1, size(markerset, 3));


nlength = squeeze(sum(all(isfinite(markerset(1:2,:,sel_markers)), 1), 2))';
if (cfg.nlength_min)
    sel_markers(nlength < cfg.nlength_minval) = false;    
end;
if (cfg.nlength_max)
    sel_markers(nlength > cfg.nlength_maxval) = false;    
end;
if (cfg.align_dist_max)
    sel_markers(aligncnf.distance_max > cfg.align_dist_maxval) = false;
end;
if (cfg.align_dist_mean)
    sel_markers(aligncnf.distance_mean > cfg.align_dist_meanval) = false;
end;
if (cfg.align_dist3d)
    sel_markers(sqrt(sum((aligncnf.X(1:3,:) - repmat(mean(aligncnf.X(1:3,:),2),[1,size(aligncnf.X,2)])) .^ 2, 1)) > cfg.align_dist3dval) = false;
end;
if (cfg.align_dist3d_regplane)
    
    X = aligncnf.X(1:3,:)';
    X = X(all(isfinite(X),2),:);

    % http://www.mathworks.com/products/statistics/demos.html?file=/products/demos/shipping/stats/orthoregdemo.html
    [coeff] = princomp(X);
    normal = coeff(:,3);

    error = abs((X - repmat(mean(X,1),size(X,1),1)) * normal)';
    
    sel_markers(error > cfg.align_dist3d_regplaneval) = false;
end;




if (~cfg.seltype_replace)
    sel_markers(handles.data.markGui_data.fcn_local.getListboxValue(handles.data.markGui_data.listbox_markers)) = true;
end;


 
handles.output.sel_markers = find(sel_markers);
handles.output.cfg = cfg;

guidata(handles.figure_markGui_SelectMarker, handles);

uiresume(handles.figure_markGui_SelectMarker);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ms_idx, ms_name, ms_fieldname, aligncnf] = getSelectedAlignment(handles);

aligncnf = [];

ss = get(handles.listbox_align_markersets, 'String');
ss2 = get(handles.listbox_align_alignments, 'String');
if (isempty(ss) || isempty(ss2))
    ms_idx = 1;
    ms_name = '';
    ms_fieldname = '';
else
    ms_idx = handles.data.align.markerset_alignfields_valid(get(handles.listbox_align_markersets, 'Value'));
    ms_name = ss{get(handles.listbox_align_markersets, 'Value')};
    ss2 = handles.data.align.markerset_alignfields{ms_idx};
    ms_fieldname = ss2{get(handles.listbox_align_alignments,'Value')};


    if (nargout >= 4)
        switch (ms_fieldname)
            case {'alignment3ditr', 'alignment3d'}
                markerset = handles.data.markGui_data.data.markersets(ms_idx);

                aligncnf.aligncnf = markerset.(ms_fieldname);
                aligncnf.markerset = handles.data.markGui_data.data.markersets(handles.data.markGui_data.data.sel_markerset).markerset;

                [aligncnf.P, xdummy, aligncnf.X, aligncnf.x] = tom_mark_cvaf_alignment3d_reproj(aligncnf.aligncnf.tiltangles, aligncnf.aligncnf.psi, aligncnf.aligncnf.tx, aligncnf.aligncnf.ty, aligncnf.aligncnf.imsize, [], aligncnf.markerset, false);
                aligncnf.x = aligncnf.x(1:2,:,:);
                aligncnf.distance = sqrt(sum((aligncnf.markerset - aligncnf.x).^2, 1));


                aligncnf.distance_mean = nan(1, size(aligncnf.distance,3));
                for (i=1:length(aligncnf.distance_mean))
                    d = aligncnf.distance(1,:,i);
                    d = d(isfinite(d));
                    if (~isempty(d))
                        aligncnf.distance_mean(i) = mean(d);
                    end;
                end;     

                aligncnf.distance_max = squeeze(max(aligncnf.distance, [], 2))';            

            otherwise
                warning('Unexpected code reached. CHECK!!!!');
        end;
    end;

end;





end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




