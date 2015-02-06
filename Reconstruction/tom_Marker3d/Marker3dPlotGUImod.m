function varargout = Marker3dPlotGUImod(varargin)
%--------------------------------------------------------------------------
% Begin initialization code - DO NOT EDIT
%--------------------------------------------------------------------------
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Marker3dPlotGUImod_OpeningFcn, ...
                   'gui_OutputFcn',  @Marker3dPlotGUImod_OutputFcn, ...
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
%--------------------------------------------------------------------------
% End initialization code - DO NOT EDIT
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% MARKER3D Plot GUI
%--------------------------------------------------------------------------
function Marker3dPlotGUImod_OpeningFcn(hObject,eventdata,handles,varargin)

handles.PlotData = varargin{1};

set(handles.dispMarkerfilename,'String',handles.PlotData.Markerfilename);
set(handles.dispMarkerfilepath,'String',handles.PlotData.Markerfilepath);
set(handles.dispV,'String',handles.PlotData.v);
set(handles.dispS,'String',handles.PlotData.s);
set(handles.dispRP,'String',handles.PlotData.rp);
set(handles.dispNoNaNs,'String',handles.PlotData.NoNaNs);
set(handles.dispImdim,'String',handles.PlotData.Imdim);

set(handles.dispRM,'String',handles.PlotData.rm);
OriginCoordinates = ...
[num2str(handles.PlotData.xrm) '  '...
 num2str(handles.PlotData.yrm) '  '...
 num2str(handles.PlotData.zrm)];
set(handles.dispOrigin,'String',OriginCoordinates);

set(handles.dispModel,'String',handles.PlotData.Model);
set(handles.dispPrecision,'String',handles.PlotData.Precision);
set(handles.dispAlgorithm,'String',handles.PlotData.Algorithm);

set(handles.dispRigidBodyDevation,'String','Not implemented.');
%set(handles.dispResiduum,'String','Not implemented.');

set(handles.dispAlphaRigidBody,'String',...
   (180/pi)*mean(handles.PlotData.alphaalg3d));

set(handles.dispMeanAlpha,'String',(180/pi)*mean(handles.PlotData.alpha));


handles.PlotData.difx = abs(handles.PlotData.tx - handles.PlotData.txalg3d);
handles.PlotData.dify = abs(handles.PlotData.ty - handles.PlotData.tyalg3d);

handles.PlotData.Meandifx = mean(handles.PlotData.difx);
handles.PlotData.Meandify = mean(handles.PlotData.dify);


set(handles.dispXdif,'String',handles.PlotData.Meandifx);
set(handles.dispYdif,'String',handles.PlotData.Meandify);




% Projection translation in x direction
axes(handles.plottx);
plot((180/pi)*handles.PlotData.theta,handles.PlotData.txalg3d,'k-');
hold on;
plot((180/pi)*handles.PlotData.theta,handles.PlotData.tx,'g-');

% Projection translation in y direction
axes(handles.plotty);
plot((180/pi)*handles.PlotData.theta,handles.PlotData.tyalg3d,'k-');
hold on;
plot((180/pi)*handles.PlotData.theta,handles.PlotData.ty,'g-');

% Tilt Axis Angle
axes(handles.plotalpha);
Meanalpha=(180/pi)*mean(handles.PlotData.alpha);
plot((180/pi)*handles.PlotData.theta,Meanalpha,'g-');
hold on;
plot((180/pi)*handles.PlotData.theta,(180/pi)*handles.PlotData.alphaalg3d,'k-');
plot((180/pi)*handles.PlotData.theta,(180/pi)*handles.PlotData.alpha,'r-');

% Projection Isoscale
axes(handles.plotisoscale);
plot((180/pi)*handles.PlotData.theta,1,'k-');
hold on;
plot((180/pi)*handles.PlotData.theta,handles.PlotData.isoscale,'b-');

% % x Scale
% axes(handles.plotxscale);
% plot((180/pi)*PlotData.theta,1,'k-');
% hold on;
% plot((180/pi)*PlotData.theta,PlotData.xscale,'b-');
% 
% % y Scale
% axes(handles.plotyscale);
% plot((180/pi)*PlotData.theta,1,'k-');
% hold on;
% plot((180/pi)*PlotData.theta,PlotData.yscale,'b-');
% 
% % z Scale
% axes(handles.plotzscale);
% plot((180/pi)*PlotData.theta,1,'k-');
% hold on;
% plot((180/pi)*PlotData.theta,PlotData.zscale,'b-');

guidata(hObject,handles);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% OUTPUT
%--------------------------------------------------------------------------
function varargout = Marker3dPlotGUImod_OutputFcn(hObject,eventdata,handles) 
%--------------------------------------------------------------------------



function plotXDif_Callback(hObject, eventdata, handles)

% if (get(hObject,'Value') == get(hObject,'Max'))
% axes(handles.plottx);
% cla;
% plot((180/pi)*handles.PlotData.theta,handles.PlotData.difx,'k-');
% 
% else
% 
% axes(handles.plottx);
% cla;
% plot((180/pi)*handles.PlotData.theta,handles.PlotData.txalg3d,'k-');
% hold on;
% plot((180/pi)*handles.PlotData.theta,handles.PlotData.tx,'g-');
% 
% end



function plotYDif_Callback(hObject, eventdata, handles)

% if (get(hObject,'Value') == get(hObject,'Max'))
% axes(handles.plotty);
% cla;
% plot((180/pi)*handles.PlotData.theta,handles.PlotData.dify,'k-');
% 
% else
% 
% axes(handles.plotty);
% cla;
% plot((180/pi)*handles.PlotData.theta,handles.PlotData.tyalg3d,'k-');
% hold on;
% plot((180/pi)*handles.PlotData.theta,handles.PlotData.ty,'g-');
% 
% end

