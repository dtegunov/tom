function tom_Marker3d(varargin)
%TOM_MARKER3D is a tool for alignment.
%
%   tom_Marker3d
%
%   You can choose between different alignment models.
%
%PARAMETERS
%
%  INPUT
%   MARKERFILE              markerfile for alignment
%
%  OUTPUT
%   ALIGNMENT              alignment parameter
%
%EXAMPLE
%
%REFERENCES
%
%SEE ALSO
%
%   created by ME 01/01/07
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



% -------------------------------------------------------------------------
% Begin initialization code - DO NOT EDIT
% -------------------------------------------------------------------------
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tom_Marker3d_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_Marker3d_OutputFcn, ...
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
% -------------------------------------------------------------------------
% End initialization code - DO NOT EDIT
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% MARKER3D Multi Model Alignment Tool
% -------------------------------------------------------------------------
function tom_Marker3d_OpeningFcn(hObject,eventdata,handles,varargin)

disp(' ');
disp('------------------------------------------------------------------');
disp('Welcome to tom_Marker3d');
disp('------------------------------------------------------------------');
disp(' ');

% Status variables
handles.Markerfileloaded = 'No';
handles.Alignmentdone = 'No';

guidata(hObject,handles);
%disp(handles);
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% OUTPUT
% -------------------------------------------------------------------------
function varargout = tom_Marker3d_OutputFcn(hObject,eventdata,handles)
% -------------------------------------------------------------------------





% *************************************************************************
% *************************************************************************
% **********     LOAD MARKERFILE     **************************************
% *************************************************************************
% *************************************************************************

% -------------------------------------------------------------------------
% Button load markerfile
% -------------------------------------------------------------------------
function LoadMarkerfile_Callback(hObject,eventdata,handles)

% Initialization
handles = Marker3dInitialisation(handles);

% Get markerfilename and path                   
[Markerfilename,Markerfilepath] = uigetfile('*.em','Locate markerfile'); 

if isequal(Markerfilename,0)|isequal(Markerfilepath,0)
    return;
else

    handles.Markerfilename = Markerfilename;
    handles.Markerfilepath = Markerfilepath;

end

% Load markerfile
workingdir = pwd;
cd (handles.Markerfilepath);
Markerfile = tom_emread(handles.Markerfilename);
cd (workingdir);

% Markerfile -> Marker3d
[theta,xm,ym,v,s,rp,NoNaNs,Imdim] = Marker3dDataCreation(Markerfile);

handles.theta = theta;
handles.xm = xm;
handles.ym = ym;
handles.v = v;
handles.s = s;
handles.rp = rp;
handles.NoNaNs = NoNaNs;
handles.Imdim = Imdim;
handles.Markerfile = Markerfile.Value;
clear Markerfilename Markerfilepath ...
      theta xm ym v s rp NoNaNs Imdim Markerfile;

% New status and display
handles.Markerfileloaded = 'Yes';
set(handles.dispMarkerfilename,'String',handles.Markerfilename);
handles.Alignmentdone = 'No'; % Reset new markerfile
set(handles.dispAlignmentdone,'String','No alignment done');

% Construct reference marker popup
for k=1:(handles.s)
    MarkerIndex(k) = k;
end
set(handles.popupReferenceMarker,'String',MarkerIndex,'Value',1);
clear MarkerIndex;

% Reset default origin
set(handles.EditXrm,'String',0); handles.xrm = 0;
set(handles.EditYrm,'String',0); handles.yrm = 0;
set(handles.EditZrm,'String',0); handles.zrm = 0;

guidata(hObject,handles);
%disp(handles);
% -------------------------------------------------------------------------




% *************************************************************************
% *************************************************************************
% **********     PANEL CHOOSE ORIGIN     **********************************
% *************************************************************************
% *************************************************************************

% -------------------------------------------------------------------------
% Popup reference marker
% -------------------------------------------------------------------------
function popupReferenceMarker_Callback(hObject,eventdata,handles)

% Choose reference marker
string_list = get(hObject,'String');
val = get(hObject,'Value');
rm = str2num(string_list(val));
handles.rm = rm;
clear string_list val rm;

% Reset default origin
set(handles.EditXrm,'String',0); handles.xrm = 0;
set(handles.EditYrm,'String',0); handles.yrm = 0;
set(handles.EditZrm,'String',0); handles.zrm = 0;

% New status
handles.Alignmentdone = 'No'; % Reset new origin
set(handles.dispAlignmentdone,'String','No alignment done');

guidata(hObject,handles);
%disp(handles);
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function popupReferenceMarker_CreateFcn(hObject,eventdata,handles)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Button default origin
% -------------------------------------------------------------------------
function DefaultOrigin_Callback(hObject,eventdata,handles)
    
% Check status
if (strcmp(handles.Markerfileloaded,'No')) == 1
    message=['No markerfile loaded.'];
    uiwait(msgbox(message,'Error','error'));
    return;
end

% Default origin
set(handles.EditXrm,'String',handles.xm(handles.rp,handles.rm));
handles.xrm = handles.xm(handles.rp,handles.rm);
set(handles.EditYrm,'String',handles.ym(handles.rp,handles.rm));
handles.yrm = handles.ym(handles.rp,handles.rm);
set(handles.EditZrm,'String',((handles.Imdim)/2+1));
handles.zrm = ((handles.Imdim)/2+1);

% New status
handles.Alignmentdone = 'No'; % Reset new origin
set(handles.dispAlignmentdone,'String','No alignment done');

guidata(hObject,handles);
%disp(handles);
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Edittext x Origin
% -------------------------------------------------------------------------
function EditXrm_Callback(hObject, eventdata, handles)

xrm = str2double(get(hObject,'string'));
if isnan(xrm)
    message=['You must enter a numeric value.'];
    uiwait(msgbox(message,'Error','error'));
end
handles.xrm = xrm;
clear xrm;

% New status
handles.Alignmentdone = 'No'; % Reset new origin
set(handles.dispAlignmentdone,'String','No alignment done');

guidata(hObject,handles);
%disp(handles);
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function EditXrm_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Edittext y Origin
% -------------------------------------------------------------------------
function EditYrm_Callback(hObject, eventdata, handles)

yrm = str2double(get(hObject,'string'));
if isnan(yrm)
    message=['You must enter a numeric value.'];
    uiwait(msgbox(message,'Error','error'));
end
handles.yrm = yrm;
clear yrm;

% New status
handles.Alignmentdone = 'No'; % Reset new origin
set(handles.dispAlignmentdone,'String','No alignment done');

guidata(hObject,handles);
%disp(handles);
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function EditYrm_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Edittext z Origin
% -------------------------------------------------------------------------
function EditZrm_Callback(hObject, eventdata, handles)

zrm = str2double(get(hObject,'string'));
if isnan(zrm)
    message=['You must enter a numeric value.'];
    uiwait(msgbox(message,'Error','error'));
end
handles.zrm = zrm;
clear zrm;

% New status
handles.Alignmentdone = 'No'; % Reset new origin
set(handles.dispAlignmentdone,'String','No alignment done');

guidata(hObject,handles);
%disp(handles);
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function EditZrm_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% *************************************************************************
% *************************************************************************
% % **********     PANEL DO ALIGNMENT     *********************************
% *************************************************************************
% *************************************************************************

% -------------------------------------------------------------------------
% Popup choose model
% -------------------------------------------------------------------------
function popupModel_Callback(hObject, eventdata, handles)

val = get(hObject,'Value');
str = get(hObject, 'String');
switch str{val};
case 'tom_alignment3d'
   handles.Model = 'tom_alignment3d';
case 'Rigid Body REF'
   handles.Model = 'Rigid Body REF';
case 'Rigid Body CM'
   handles.Model = 'Rigid Body CM';
case 'Free Tilt REF'
   handles.Model = 'Free Tilt REF';
case 'Free Tilt CM'
   handles.Model = 'Free Tilt CM';
case 'Free Tilt IsoScale REF'
   handles.Model = 'Free Tilt IsoScale REF';
case 'Free Tilt IsoScale CM'
   handles.Model = 'Free Tilt IsoScale CM';
case 'Free Tilt 3dScale REF'
   handles.Model = 'Free Tilt 3dScale REF';
case 'Free Tilt 3dScale CM'
   handles.Model = 'Free Tilt 3dScale CM';
end

% New status
handles.Alignmentdone = 'No'; % Reset new model
set(handles.dispAlignmentdone,'String','No alignment done');

guidata(hObject,handles);
%disp(handles);
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function popupModel_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Popup choose precision
% -------------------------------------------------------------------------
function popupPrecision_Callback(hObject, eventdata, handles)

val = get(hObject,'Value');
str = get(hObject, 'String');
switch str{val};
case '0.001'
   handles.Precision = '0.001';
case '0.01'
   handles.Precision = '0.01';
case '0.1'
   handles.Precision = '0.1';
end

% New status
handles.Alignmentdone = 'No'; % Reset new precision
set(handles.dispAlignmentdone,'String','No alignment done');

guidata(hObject,handles);
%disp(handles);
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function popupPrecision_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Popup choose algorithm
% -------------------------------------------------------------------------
function popupAlgorithm_Callback(hObject, eventdata, handles)

val = get(hObject,'Value');
str = get(hObject, 'String');
switch str{val};
case 'Conjugate gradients'
   handles.Algorithm = 'Conjugate gradients';
case 'Simplex search'
   handles.Algorithm = 'Simplex search';
end

% New status
handles.Alignmentdone = 'No'; % Reset new algorithm
set(handles.dispAlignmentdone,'String','No alignment done');

guidata(hObject,handles);
%disp(handles);
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function popupAlgorithm_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Button do alignment now
% -------------------------------------------------------------------------
function DoAlignmentNow_Callback(hObject,eventdata,handles,v,s,xm,ym,theta)

% Check status
if (strcmp(handles.Markerfileloaded,'No')) == 1
    message=['No markerfile loaded.'];
    uiwait(msgbox(message,'Error','error'));
    return;
end

% ********************************************************************
% **********     Start values     ****************************************
% ********************************************************************
[Matrixmark,psi,sigma,x,y,z] = tom_alignment3d...
(handles.Markerfile,handles.rm,handles.rp,...
[handles.xrm handles.yrm handles.zrm],handles.Imdim);

handles.x3dalg3d = x; 
handles.y3dalg3d = y;
handles.z3dalg3d = z;
handles.txalg3d = Matrixmark(7,:,1);
handles.tyalg3d = Matrixmark(8,:,1);
handles.alphaalg3d(1:handles.v) = psi;

clear Matrixmark psi sigma x y z;
% ********************************************************************

% ********************************************************************
% **********     tom_alignment3d alignment model     ***************
% ********************************************************************
if (strcmp(handles.Model,'tom_alignment3d')) == 1 

% copy start values to handles
handles.Residuum = 0;

handles.x3d = handles.x3dalg3d+(handles.Imdim/2+1)-(handles.xrm);
handles.y3d = handles.y3dalg3d+(handles.Imdim/2+1)-(handles.yrm);
handles.z3d = handles.z3dalg3d+(handles.Imdim/2+1)-(handles.zrm);

handles.tx = handles.txalg3d;
handles.ty = handles.tyalg3d;
handles.alpha = handles.alphaalg3d;
handles.isoscale(1:handles.v) = 1;
handles.xscale(1:handles.v) = 1;
handles.yscale(1:handles.v) = 1;
handles.zscale(1:handles.v) = 1;

% alignment done
message=['tom_alignment3d alignment done .'];
uiwait(msgbox(message,'Status','warn'));

end
% ********************************************************************

% ********************************************************************
% **********     Rigid Body REF alignment model     ******************
% ********************************************************************
if (strcmp(handles.Model,'Rigid Body REF')) == 1  

% optimization
[Residuum,x3d,y3d,z3d,tx,ty,alpha,...
 isoscale,xscale,yscale,zscale] = Marker3dRigidBodyREF(handles);

% copy results to handles
handles.Residuum = Residuum;
handles.x3d = x3d;
handles.y3d = y3d;
handles.z3d = z3d;
handles.tx = tx;
handles.ty = ty;
handles.alpha = alpha;
handles.isoscale = isoscale;
handles.xscale = xscale;
handles.yscale = yscale;
handles.zscale = zscale;

clear Residuum x3d y3d z3d tx ty alpha isoscale xscale yscale zscale;

% Message Rigid Body REF alignment done
message=['Rigid Body REF alignment done .'];
uiwait(msgbox(message,'Status','warn'));

end
% ********************************************************************

% ********************************************************************
% **********     Rigid Body CM alignment model     *******************
% ********************************************************************
if (strcmp(handles.Model,'Rigid Body CM')) == 1  

% optimization
[Residuum,x3d,y3d,z3d,tx,ty,alpha,...
 isoscale,xscale,yscale,zscale] = Marker3dRigidBodyCM(handles);


% copy results to handles
handles.Residuum = Residuum;
handles.x3d = x3d;
handles.y3d = y3d;
handles.z3d = z3d;
handles.tx = tx;
handles.ty = ty;
handles.alpha = alpha;
handles.isoscale = isoscale;
handles.xscale = xscale;
handles.yscale = yscale;
handles.zscale = zscale;

clear Residuum x3d y3d z3d tx ty alpha isoscale xscale yscale zscale;

% Message Rigid Body CM alignment done
message=['Rigid Body CM alignment done .'];
uiwait(msgbox(message,'Status','warn'));

end
% ********************************************************************

% ********************************************************************
% ***********     Free Tilt REF alignment model     ********************
% ********************************************************************
if (strcmp(handles.Model,'Free Tilt REF')) == 1  

% optimization
[Residuum,x3d,y3d,z3d,tx,ty,alpha,...
             isoscale,xscale,yscale,zscale] = Marker3dFreeTiltREF(handles);


% copy results to handles
handles.Residuum = Residuum;
handles.x3d = x3d;
handles.y3d = y3d;
handles.z3d = z3d;
handles.tx = tx;
handles.ty = ty;
handles.alpha = alpha;
handles.isoscale = isoscale;
handles.xscale = xscale;
handles.yscale = yscale;
handles.zscale = zscale;

clear Residuum x3d y3d z3d tx ty alpha isoscale xscale yscale zscale;

% Message Free Tilt REF alignment done
message=['Free Tilt REF alignment done .'];
uiwait(msgbox(message,'Status','warn'));

end
% ********************************************************************

% ********************************************************************
% ***********     Free Tilt CM alignment model     ********************
% ********************************************************************
if (strcmp(handles.Model,'Free Tilt CM')) == 1  

% optimization
[Residuum,x3d,y3d,z3d,tx,ty,alpha,...
              isoscale,xscale,yscale,zscale] = Marker3dFreeTiltCM(handles);


% copy results to handles
handles.Residuum = Residuum;
handles.x3d = x3d;
handles.y3d = y3d;
handles.z3d = z3d;
handles.tx = tx;
handles.ty = ty;
handles.alpha = alpha;
handles.isoscale = isoscale;
handles.xscale = xscale;
handles.yscale = yscale;
handles.zscale = zscale;

clear Residuum x3d y3d z3d tx ty alpha isoscale xscale yscale zscale;

% Message Free Tilt CM alignment done
message=['Free Tilt CM alignment done .'];
uiwait(msgbox(message,'Status','warn'));

end
% ********************************************************************

% ********************************************************************
% ********     Free Tilt IsoScale REF alignment model     *****************
% ********************************************************************
if (strcmp(handles.Model,'Free Tilt IsoScale REF')) == 1  

% optimization
[Residuum,x3d,y3d,z3d,tx,ty,alpha,...
     isoscale,xscale,yscale,zscale] = Marker3dFreeTiltIsoScaleREF(handles);

 
% copy results to handles
handles.Residuum = Residuum;
handles.x3d = x3d;
handles.y3d = y3d;
handles.z3d = z3d;
handles.tx = tx;
handles.ty = ty;
handles.alpha = alpha;
handles.isoscale = isoscale;
handles.xscale = xscale;
handles.yscale = yscale;
handles.zscale = zscale;

clear Residuum x3d y3d z3d tx ty alpha isoscale xscale yscale zscale;

% Message Free Tilt IsoScale REF alignment done
message=['Free Tilt IsoScale REF alignment done .'];
uiwait(msgbox(message,'Status','warn'));

end
% ********************************************************************

% ********************************************************************
% ********     Free Tilt IsoScale CM alignment model     *****************
% ********************************************************************
if (strcmp(handles.Model,'Free Tilt IsoScale CM')) == 1  

% optimization
[Residuum,x3d,y3d,z3d,tx,ty,alpha,...
       isoscale,xscale,yscale,zscale] = Marker3dFreeTiltIsoScaleCM(handles);


% copy results to handles
handles.Residuum = Residuum;
handles.x3d = x3d;
handles.y3d = y3d;
handles.z3d = z3d;
handles.tx = tx;
handles.ty = ty;
handles.alpha = alpha;
handles.isoscale = isoscale;
handles.xscale = xscale;
handles.yscale = yscale;
handles.zscale = zscale;

clear Residuum x3d y3d z3d tx ty alpha isoscale xscale yscale zscale;

% Message Free Tilt IsoScale CM alignment done
message=['Free Tilt IsoScale CM alignment done .'];
uiwait(msgbox(message,'Status','warn'));

end
% ********************************************************************

% ********************************************************************
% ********     Free Tilt 3dScale REF alignment model     *****************
% ********************************************************************
if (strcmp(handles.Model,'Free Tilt 3dScale REF')) == 1  

% optimization
[Residuum,x3d,y3d,z3d,tx,ty,alpha,...
      isoscale,xscale,yscale,zscale] = Marker3dFreeTilt3dScaleREF(handles);


% copy results to handles
handles.Residuum = Residuum;
handles.x3d = x3d;
handles.y3d = y3d;
handles.z3d = z3d;
handles.tx = tx;
handles.ty = ty;
handles.alpha = alpha;
handles.isoscale = isoscale;
handles.xscale = xscale;
handles.yscale = yscale;
handles.zscale = zscale;

clear Residuum x3d y3d z3d tx ty alpha isoscale xscale yscale zscale;

% Message Free Tilt 3dScale REF alignment done
message=['Free Tilt 3dScale REF alignment done .'];
uiwait(msgbox(message,'Status','warn'));

end
% ********************************************************************

% ********************************************************************
% ********     Free Tilt 3dScale CM alignment model     *****************
% ********************************************************************
if (strcmp(handles.Model,'Free Tilt 3dScale CM')) == 1  

% optimization
[Residuum,x3d,y3d,z3d,tx,ty,alpha,...
        isoscale,xscale,yscale,zscale] = Marker3dFreeTilt3dScaleCM(handles);


% copy results to handles
handles.Residuum = Residuum;
handles.x3d = x3d;
handles.y3d = y3d;
handles.z3d = z3d;
handles.tx = tx;
handles.ty = ty;
handles.alpha = alpha;
handles.isoscale = isoscale;
handles.xscale = xscale;
handles.yscale = yscale;
handles.zscale = zscale;

clear Residuum x3d y3d z3d tx ty alpha isoscale xscale yscale zscale;

% Message Free Tilt 3dScale CM alignment done
message=['Free Tilt 3dScale CM alignment done .'];
uiwait(msgbox(message,'Status','warn'));

end
% ********************************************************************


% New status
handles.Alignmentdone = 'Yes';
set(handles.dispAlignmentdone,'String',handles.Markerfilename);

guidata(hObject,handles);
%disp(handles);
% -------------------------------------------------------------------------


% *************************************************************************
% *************************************************************************
% PANEL PLOT AND SAVE ALIGNMENT
% *************************************************************************
% *************************************************************************

% -------------------------------------------------------------------------
% Popup tiltseries
% -------------------------------------------------------------------------
function popupTiltSeries_Callback(hObject,eventdata,handles)

% Choose tiltseries
val = get(hObject,'Value');
str = get(hObject, 'String');
switch str{val};
case '0'
   handles.Tiltseries = '0';
case '1'
   handles.Tiltseries = '1';
case '2'
   handles.Tiltseries = '2';
end

guidata(hObject,handles);
%disp(handles);
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function popupTiltSeries_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Button plot alignment
% -------------------------------------------------------------------------
function PlotAlignment_Callback(hObject,eventdata,handles)

% Check status
if (strcmp(handles.Markerfileloaded,'No')) == 1
    message=['No alignment done.'];
    uiwait(msgbox(message,'Error','error'));
    return;
end

if (strcmp(handles.Alignmentdone,'No')) == 1
    message=['No alignment done.'];
    uiwait(msgbox(message,'Error','error'));
    return;
end

% Plot alignment GUI
inputPlotGUI = handles;
Marker3dPlotGUImod(inputPlotGUI);
clear inputPlotGUI;
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Button save alignment
% -------------------------------------------------------------------------
function SaveAlignment_Callback(hObject, eventdata, handles)

% Check status
if (strcmp(handles.Markerfileloaded,'No')) == 1
    message=['No alignment done.'];
    uiwait(msgbox(message,'Error','error'));
    return;
end

if (strcmp(handles.Alignmentdone,'No')) == 1
    message=['No alignment done.'];
    uiwait(msgbox(message,'Error','error'));
    return;
end

% Save alignment of tiltserie 0 to workspace
if (strcmp(handles.Tiltseries,'0')) == 1

    assignin('base','Tiltangle',(180/pi)*handles.theta);

    handles.m3d(1,1:handles.s) = handles.x3d(1:handles.s);
    handles.m3d(2,1:handles.s) = handles.y3d(1:handles.s);
    handles.m3d(3,1:handles.s) = handles.z3d(1:handles.s);
    assignin('base','m3d',handles.m3d);

    % Substract Imdim
    tx = handles.tx - handles.Imdim./2;
    assignin('base','tx',tx);
    
    ty = handles.ty - handles.Imdim./2;
    assignin('base','ty',ty);
    
    % Add 90 degrees
    alpha = (180/pi)*handles.alpha + 90;
    assignin('base','Tiltaxis',alpha);
    
    assignin('base','isoscale',handles.isoscale);

end

% Save alignment of tiltserie 1 to workspace
if (strcmp(handles.Tiltseries,'1')) == 1

    assignin('base','Tiltangle1',(180/pi)*handles.theta);

    handles.m3d(1,1:handles.s) = handles.x3d(1:handles.s);
    handles.m3d(2,1:handles.s) = handles.y3d(1:handles.s);
    handles.m3d(3,1:handles.s) = handles.z3d(1:handles.s);
    assignin('base','m3d1',handles.m3d);

    assignin('base','tx1',handles.tx);
    assignin('base','ty1',handles.ty);
    assignin('base','Tiltaxis1',(180/pi)*handles.alpha);

end

% Save alignment of tiltserie 2 to workspace
if (strcmp(handles.Tiltseries,'2')) == 1

    assignin('base','Tiltangle2',(180/pi)*handles.theta);

    handles.m3d(1,1:handles.s) = handles.x3d(1:handles.s);
    handles.m3d(2,1:handles.s) = handles.y3d(1:handles.s);
    handles.m3d(3,1:handles.s) = handles.z3d(1:handles.s);
    assignin('base','m3d2',handles.m3d);

    assignin('base','tx2',handles.tx);
    assignin('base','ty2',handles.ty);
    assignin('base','Tiltaxis2',(180/pi)*handles.alpha);

end

message=['Alignment saved to workspace.'];
uiwait(msgbox(message,'Message','warn'));
return;
% -------------------------------------------------------------------------



% *************************************************************************
% *************************************************************************
% **********     INCLUDED SUBFUNCTIONS     ********************************
% *************************************************************************
% *************************************************************************

% -------------------------------------------------------------------------
% INITIALIZATION
% -------------------------------------------------------------------------
function handles = Marker3dInitialisation(handles)

handles.Markerfilename = '';
handles.Markerfilepath = '';
handles.Alignmentfilename = '';
handles.Alignmentfilepath = '';

handles.theta = 0;        % Measured tilt angles
handles.xm = 0;           % Measured x coordinates
handles.ym = 0;           % Measured y coordinates 
handles.v = 0;            % Number of projections
handles.s = 0;            % Number of markers on each projection
handles.rp = 0;           % Minimum tilt projection (zero degree)
handles.NoNaNs = 0;       % Number of data points not measured
handles.Imdim = 0;        % Image dimension
handles.Markerfile = 0;   % Markerfile Value for start values and plot

handles.rm = 1;           % Default reference marker
handles.xrm = 0;          % Default x coordinates of reference marker
handles.yrm = 0;          % Default y coordinates of reference marker
handles.zrm = 0;          % Default z coordinates of reference marker

handles.Model = 'tom_alignment3d'; % Which optimization model
handles.Precision = '0.001'; % Which optimization precision
handles.Algorithm = 'Conjugate gradients'; % Which optimization algorithm

handles.Residuum = 0;     % RESIDUUM
handles.x3d = 0;          % x values of 3d marker coordinates
handles.y3d = 0;          % y values of 3d marker coordinates
handles.z3d = 0;          % z values of 3d marker coordinates
handles.m3d = 0;          % x y z values combined

handles.tx = 0;           % Projection x translation
handles.ty = 0;           % Projection y translation
handles.alpha = 0;        % Tilt axis angle
handles.isoscale = 1;     % Projection isoscale
handles.xscale = 1;       % Scale x
handles.yscale = 1;       % Scale y
handles.zscale = 1;       % Scale z

handles.Tiltseries = '0'; % Which tiltseries

handles.x3dalg3d = 0;     % tom_alignment3d x3d / Start values 
handles.y3dalg3d = 0;     % tom_alignment3d y3d / Start values 
handles.z3dalg3d = 0;     % tom_alignment3d z3d / Start values 
handles.txalg3d = 0;      % tom_alignment3d x translation / Start values
handles.tyalg3d = 0;      % tom_alignment3d y translation / Start values
handles.alphaalg3d = 0;   % tom_alignment3d Tilt axis angle / Start values
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% DATA CREATION
%--------------------------------------------------------------------------
function [theta,xm,ym,v,s,rp,NoNaNs,Imdim] = Marker3dDataCreation(Markerfile);

% Dimensions markerfile
v = size(Markerfile.Value,2);  % Number of projections
s = size(Markerfile.Value,3);  % Number of markers on each projection

% Real data input
xm(1:v,1:s) = Markerfile.Value(2,1:v,1:s);  % Measured x coordinates
ym(1:v,1:s) = Markerfile.Value(3,1:v,1:s);  % Measured y coordinates 
theta(1:v) = (pi./180).*Markerfile.Value(1,1:v,1);  % Measured tilt angles

% Minimum tilt projection (reference projection)
[notneeded,rp] = min(abs(theta));

% Image dimension
if size(Markerfile.Value,1) == 3
    Imdim = 1024; % Old markerfiles only have 10 lines.
else
Imdim = Markerfile.Value(12,1,1);
end

% Change (-1) into NaN (Not measured data point)
for j=1:s
   for i=1:v
      if xm(i,j) == (-1)
         xm(i,j) = NaN;
         ym(i,j) = NaN;
      end
   end
end

% Count not measured data point
zaehler=0;
for j=1:s
   for i=1:v
      if isnan(xm(i,j)) == 1
         zaehler=zaehler+1;
      end
   end
end
NoNaNs=2*zaehler;
%--------------------------------------------------------------------------


