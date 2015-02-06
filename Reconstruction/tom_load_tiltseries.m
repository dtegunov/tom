function varargout = tom_load_tiltseries(varargin)
%TOM_LOAD_TILTSERIES can just be used by TOM_FIG_MENU (TOM_SETMARK)
%
%   varargout = tom_load_tiltseries(varargin)
%
%   This routine open a dialog box window to help the user to open a
%   tilt serie. There is several field to fill manually. The user can
%   also use the 'Browse' button to fill all the field automatically. One
%   just has to select the last picture of his tilt serie.
%   ! Do not use this function if you don't use TOM_FIG_MENU
%
%PARAMETERS
%
%  INPUT
%   varargin            ...
%  
%  OUTPUT
%   varargout   		...
%
%EXAMPLE
%   ... = tom_amira_createisosurface(...);
%   creates ...
%
%REFERENCES
%
%SEE ALSO
%   TOM_FIG_MENU
%
%   created by WDN 08/12/02
%   updated by WDN 02/17/04
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

if nargin == 0  % LAUNCH GUI

	fig = openfig(mfilename,'reuse');
    
	% Use system color scheme for figure:
	set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));
	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
	guidata(fig, handles)
	if nargout > 0
		varargout{1} = fig;
    end

elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

	try
		if (nargout)
			[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
        else
			feval(varargin{:}); % FEVAL switchyard
        end
    catch
		disp(lasterr);
    end

end




% --------------------------------------------------------------------
function varargout = lb_ok_Callback(h, eventdata, handles, varargin)

param.mypathname=get(handles.lb_path,'String');
param.myfilename=get(handles.lb_file,'String');
param.myfirstnb=get(handles.lb_firstnb,'String');
param.mylastnb=get(handles.lb_lastnb,'String');
param.myext=get(handles.lb_ext,'String');
if ~isempty(param.myext)
    param.myext=['.' param.myext];
end
param.myfilemarker_default=[param.mypathname get(handles.lb_filemarker,'String')];
param.myfilemarker='';
param.Matrixmark=-1*ones(12,str2num(param.mylastnb),1);%define matrix(12,number of image; number of marks;)
param.Matrixmark(4:12,:,1)= 0;%(4:10,:,1) %ligne 1:= tilt angle ligne 2:= x of the mark ligne 3:= y of the mark
aa=[param.mypathname param.myfilename param.myfirstnb param.myext];
ah=tom_reademheader(aa);
tilt(1)=ah.Header.Parameter(19)./1000.0;
param.Matrixmark(12,1,1)=ah.Header.Size(1);
aa=[param.mypathname param.myfilename '2' param.myext];
ah=tom_reademheader(aa);
tilt(2)=ah.Header.Parameter(19)./1000.0;
param.Matrixmark(12,2,1)=ah.Header.Size(1);
if abs(tilt(1)) < abs(tilt(2))
    little=abs(tilt(1));
    No_ima=1;
else
    little=abs(tilt(2));
    No_ima=2;
end
for i=3:str2num(param.mylastnb)                        
    aa=[param.mypathname param.myfilename num2str(i) param.myext];
    ah=tom_reademheader(aa);
    param.Matrixmark(12,i,1)=ah.Header.Size(1);
    tilt(i)=ah.Header.Parameter(19)./1000.0;
    if abs(tilt(i)) < little
        little=abs(tilt(i));
        No_ima=i;
    end
end
param.angle=tilt;
param.image_ref=num2str(No_ima);
%param.Matrixmark=-1*ones(10,str2num(param.mylastnb),1);%define matrix(10,number of image; number of marks;)
                                                        %ligne 1:= tilt angle ligne 2:= x of the mark ligne 3:= y of the mark
param.Matrixmark(1,:,:)=tilt;
for i=str2num(param.myfirstnb):size(param.Matrixmark,2)
    param.Matrixmark(11,i,:)=i;
end
param.newproj_cancel='no';
param.Projection_size=ah.Header.Size(1);
%case of tom_setmark call tom_load
hh=findobj('Tag','menu');
set(hh,'Userdata',param);
close;
% --------------------------------------------------------------------
function varargout = lb_browse_Callback(h, eventdata, handles, varargin)
set(handles.lb_path,'String','');
set(handles.lb_file,'String','');
set(handles.lb_firstnb,'String','');
set(handles.lb_lastnb,'String','');
set(handles.lb_ext,'String','');
set(handles.lb_filemarker,'String','');

[f, p] = uigetfile('*.*', '---- Click on the last number of your tilt series ----');
if isequal(f,0) | isequal(p,0) 
    %error('No data loaded.'); 
    return;%nothing because Cancel is ckicked
else
    myext='';myfile='';mynb='';        
    p=strrep(p,'\','/');
    set(handles.lb_path,'String',p);
    fext=findstr(f,'.');
    if ~isempty(fext)
        fext=fext(size(fext,2));
        for i=fext+1:size(f,2)
            myext=strcat(myext,f(i));
        end
    else
        fext=size(f,2)+1;
    end
    set(handles.lb_ext,'String',myext);
    fnb=findstr(f,'_');
    if ~isempty(fnb)
        a=size(fnb,2);
        for i=fnb(a)+1:fext-1
            mynb=strcat(mynb,f(i));
        end
    end
    set(handles.lb_lastnb,'String',mynb);
    for i=1:fnb(a)
        myfile=strcat(myfile,f(i));
    end
    set(handles.lb_file,'String',myfile);
    set(handles.lb_firstnb,'String','1');
    set(handles.lb_filemarker,'String','markfile_temp.em');
end
% --------------------------------------------------------------------
function varargout = lb_cancel_Callback(h, eventdata, handles, varargin)
hh=findobj('Tag','menu');
param=get(hh,'Userdata');
param.newproj_cancel='yes';
hh=findobj('Tag','menu');
set(hh,'Userdata',param);
close;

