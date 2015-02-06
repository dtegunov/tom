function varargout = tom_os3_filterSetup(varargin)
% TOM_OS3_FILTERSETUP M-file for tom_os3_filterSetup.fig
%      TOM_OS3_FILTERSETUP, by itself, creates a new TOM_OS3_FILTERSETUP or raises the existing
%      singleton*.
%
%      H = TOM_OS3_FILTERSETUP returns the handle to a new TOM_OS3_FILTERSETUP or the handle to
%      the existing singleton*.
%
%      TOM_OS3_FILTERSETUP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TOM_OS3_FILTERSETUP.M with the given input arguments.
%
%      TOM_OS3_FILTERSETUP('Property','Value',...) creates a new TOM_OS3_FILTERSETUP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before tom_os3_display_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to tom_os3_filterSetup_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help tom_os3_filterSetup

% Last Modified by GUIDE v2.5 14-Mar-2008 14:38:46

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tom_os3_filterSetup_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_os3_filterSetup_OutputFcn, ...
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
% End initialization code - DO NOT EDIT


% --- Executes just before tom_os3_filterSetup is made visible.
function tom_os3_filterSetup_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to tom_os3_filterSetup (see VARARGIN)

% Choose default command line output for tom_os3_filterSetup
handles.output = hObject;

set(handles.rbtn_image,'Value',1);
set(handles.rbtn_peaks,'Value',0);
% Update handles structure

handles.filter.combinations = {};
handles.options = [];
guidata(hObject, handles);

% UIWAIT makes tom_os3_filterSetup wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = tom_os3_filterSetup_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function axisButtonDownCallback(hObject, eventdata, handles)

position = get(get(hObject,'Parent'),'CurrentPoint');

x = ceil(position(1));
y = ceil(position(3));

imageAxis    = handles.imageAxis;
axes(imageAxis);
hold on;
plot(x,y,'o','Color',[1 0 0]);
text(x+5,y+5,num2str(handles.img.counter),'Color',[1 0.1 0.1]);
axis off;

% disp([ num2str(handles.img.counter) ' - CCC: ' num2str(handles.img.result.ccc(x,y)) ' FING: ' num2str(handles.img.result.angleCorr(x,y)) ' Prod:' num2str(handles.img.result.autoc(x,y) * handles.img.result.angleCorr(x,y)) ]);
handles.img.counter = handles.img.counter +1;

if(isfield(handles.img,'picks'))
    picks = handles.img.picks;
    picks{length(picks)+1} = [x y];
else
    picks = {[x y]};
end;


handles.img.picks = picks;

guidata(hObject, handles);

% --------------------------------------------------------------------
function load_Callback(hObject, eventdata, handles)
% hObject    handle to load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(isfield(handles,'img'))
    if(isfield(handles.img,'training'))
        t = handles.img.training;
    end;
    handles= rmfield(handles,'img');
    if(exist('t'))
        handles.img.training = t;
    end;
end;
if(isfield(handles,'path'))
    path = handles.img.path;
else
    path = '.';
end;
[FileName,PathName,FilterIndex] = uigetfile('.mat',path);
handles.img.path = PathName;
load([PathName FileName]);

%%
peaks = get_peaks(result);


%%
result.job.volumeFile
volume = tom_emreadc3(result.job.volumeFile);
volume = volume.Value;

if(isnan(result.job.options.modifications.binning))
    result.job.options.modifications.binning = 0;
end;



%%
imageAxis    = handles.imageAxis;
colormap gray;
axes(imageAxis);
volume=tom_bin(volume,result.job.options.modifications.binning);
volume = tom_filter(volume,2);
imagesc(volume');
drawnow;

%% filter image?
filter = inputdlg('Filter image?');
if(numel(filter) > 0)
    try
       volume = tom_filter(volume,str2double(filter{1}));
       imagesc(volume');
    catch
    end;
end;
%%
set(get(imageAxis, 'Children'),'ButtonDownFcn','tom_os3_filterSetup(''axisButtonDownCallback'',gcbo,[],guidata(gcbo))')
axis off;


%%
handles.img.peaks = peaks;
if(isnan(result.job.options.modifications.binning))
    result.job.options.modification.binning = 1;
end;
handles.img.result = result;
handles.img.counter = 1;
handles.img.volume = volume;
handles.img.filename = [PathName FileName];
if(isfield(result,'angCorr'))
    handles.img.angCorrelation = result.angCorr;
end;
guidata(hObject, handles);




% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_3_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)







% --------------------------------------------------------------------
function adjustImageContrast_Callback(hObject, eventdata, handles)
% hObject    handle to adjustImageContrast (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
imageAxis=handles.imageAxis;
imcontrast(imageAxis);



% --------------------------------------------------------------------
function reload_Callback(hObject, eventdata, handles)
% hObject    handle to reload (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

load(handles.img.filename);

peaks = get_peaks(result);

volume = tom_emread(result.job.volumeFile);
volume = volume.Value;

if(isnan(result.job.options.modifications.binning))
    result.job.options.modifications.binning = 0;
end;

imageAxis    = handles.imageAxis;
colormap gray;
axes(imageAxis);
volume=tom_bin(volume,result.job.options.modifications.binning);
imagesc(volume');

%% filter image?
filter = inputdlg('Filter image?');
if(numel(filter) > 0)
    try
       volume = tom_filter(volume,str2double(filter{1}));
       imagesc(volume');
    catch
    end;
end;
%%
set(get(imageAxis, 'Children'),'ButtonDownFcn','tom_os3_display(''axisButtonDownCallback'',gcbo,[],guidata(gcbo))')
axis off;

handles.img.peaks = peaks;
handles.img.result = result;
handles.img.counter = 1;
handles.img.volume = volume;
if(isfield(result,'angCorr'))
    handles.img.angCorrelation = result.angCorr;
end;
guidata(hObject, handles);




% --------------------------------------------------------------------
function classifyStack_Callback(hObject, eventdata, handles)
% hObject    handle to classifyStack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


    if(~isfield(handles.img,'training'))
        
        %if no training values are given
        numberOfEigenImages=str2double(INPUTDLG('How many Eigenimages?'));
        numberOfCentroidImages=str2double(INPUTDLG('How many Centroid Images?'));    
        mask = tom_os3_sphereMask(particleStack(:,:,1));
        training = [];
        alignment = [];    
        normflag = 1;
    else
        
        %if training is given, use these values
        numberOfEigenImages = handles.img.training.numberEigenImages;
        numberOfCentroidImages = 0;
        mask = handles.img.training.structure.template(:,:,2);
        training = handles.img.training.values;
        alignment = handles.img.training.structure.alignment;
        normflag = 0;
    end;
    
    
    particleStack = handles.img.particleStack.stack;
    
    %normalize the stack under the mask
    mask = tom_os3_sphereMask(zeros([size(particleStack,1),size(particleStack,2)]));
    for i=1:size(particleStack,3)
        p(:,:,i)= tom_bin(tom_norm(particleStack(:,:,i),'phase',mask),1);
    end;
    particleStack = p;
    


    sigma =str2double(INPUTDLG('Maximal cluster distance in STDs. 0 - 5 is reasonable.'));    

    classifiedParticleStack = tom_os3_classifyImageStack(particleStack,handles.img.result.job.dimension,numberOfEigenImages,numberOfCentroidImages,mask,normflag,sigma,training,alignment);


projAxis = handles.projectionAxis;
axes(projAxis);
sumStack = sum(classifiedParticleStack,3);
imagesc(sumStack .* (sumStack > -10));
% axes(handles.eigenAxis1);
% imagesc(reshape(coeff(1,:),size(classfiedParticleStack,1),size(classfiedParticleStack,2)));
figure;tom_dspcub(classifiedParticleStack);

handles.img.classifiedParticleStack = classifiedParticleStack;
guidata(hObject, handles);


% --------------------------------------------------------------------
function loadTrainingStack_Callback(hObject, eventdata, handles)
% hObject    handle to loadTrainingStack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


  
    [FileName,PathName,FilterIndex] = uigetfile('.mat',path);
    
    load([PathName  FileName]);
    
%%    
    numberEigenImages=str2double(INPUTDLG('How many EigenImages?'));
    
    [a b trainingValues]=tom_os3_classifyImageStack(training.stack,2,numberEigenImages,size(training.stack,3),training.template(:,:,2),0,0,[],training.alignment);
    
    handles.img.training.values= trainingValues;
    handles.img.training.structure = training;
    handles.img.training.numberEigenImages = numberEigenImages;
    guidata(hObject, handles);

% --------------------------------------------------------------------   
function   peaks= get_peaks(result)

% peaks = result.ccc + 2*result.angleCorr;
% peaks = result.psr;
% peaks = result.autoc;
% peaks = result.angleCorr;
% peaks = result.psr .* result.autoc .* result.angleCorr;
% peaks = result.autoc .* result.angleCorr;
% peaks = result.ccc .* result.psr .* result.autoc .* result.angleCorr;
%   peaks = result.psr .* result.angleCorr;
% peaks = result.ccc .* result.psr .* result.autoc;

%7mu best combination
a = 0.55; %0.89
b = 0.11; %0.07
c = 0.9;  %0.51
peaks = a * result.ccc + b * result.psr + c * result.autoc;


%20S combination

% a = 0.58; %0.89
% b = 0.34; %0.07
% c = 0.17;  %0.51
% peaks = a * result.ccc + b * result.psr + c * result.autoc;


% %contest combination
% a = 0.75; %0.89
% b = 0.84; %0.07
% c = 0.32;  %0.51
% peaks = a * result.ccc + b * result.psr + c * result.autoc;
% 
% %2712
a = 0.21; %0.89
b = 0.05; %0.07
c = 0.94;  %0.51
peaks = a * result.ccc + b * result.psr + c * result.autoc;


%phantom
a = 0; %0.89
b = 0; %0.07
c = 1;  %0.51
peaks = a * result.ccc + b * result.psr + c * result.autoc;
% --------------------------------------------------------------------
function Untitled_19_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Load_Picklist_Callback(hObject, eventdata, handles)
% hObject    handle to Load_Picklist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
result = handles.img.result;

[FileName,PathName,FilterIndex] = uigetfile('.mat',path);
handles.img.path = PathName;
load([PathName FileName]);


number=str2double(INPUTDLG('How many?'));

for i=1:number
    
    handPick = align2d(i);
    x = handPick.position.x;
    y = handPick.position.y;
    
    
    axes(handles.imageAxis);
    hold on;
    plot(x,y,'o','Color',[1 0 0],'MarkerSize',result.job.templateSize(1)/2^(result.job.options.modifications.binning+1),'LineWidth',2);
    text(x+5,y+5,num2str(i),'Color',[1 0 0]);
    hold off;
    axis off;
    
end;


% --- Executes on slider movement.
function xcfSlider_Callback(hObject, eventdata, handles)
% hObject    handle to xcfSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
sliderValue = get(hObject,'Value');

ed = handles.xcfEdit;

set(ed,'String',num2str(sliderValue));
handles.img.coeffs.xcf = sliderValue;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function xcfSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xcfSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function socSlider_Callback(hObject, eventdata, handles)
% hObject    handle to socSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
sliderValue = get(hObject,'Value');

ed = handles.socEdit;

set(ed,'String',num2str(sliderValue));
handles.img.coeffs.soc = sliderValue;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function socSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to socSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function psrSlider_Callback(hObject, eventdata, handles)
% hObject    handle to psrSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
sliderValue = get(hObject,'Value');

ed = handles.psrEdit;

set(ed,'String',num2str(sliderValue));
handles.img.coeffs.psr = sliderValue;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function psrSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to psrSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function xcfEdit_Callback(hObject, eventdata, handles)
% hObject    handle to xcfEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xcfEdit as text
%        str2double(get(hObject,'String')) returns contents of xcfEdit as a double


% --- Executes during object creation, after setting all properties.
function xcfEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xcfEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function psrEdit_Callback(hObject, eventdata, handles)
% hObject    handle to psrEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of psrEdit as text
%        str2double(get(hObject,'String')) returns contents of psrEdit as a double


% --- Executes during object creation, after setting all properties.
function psrEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to psrEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function socEdit_Callback(hObject, eventdata, handles)
% hObject    handle to socEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of socEdit as text
%        str2double(get(hObject,'String')) returns contents of socEdit as a double


% --- Executes during object creation, after setting all properties.
function socEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to socEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btnPick.
function btnPick_Callback(hObject, eventdata, handles)
% hObject    handle to btnPick (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
result = handles.img.result;

xcf = handles.img.coeffs.xcf;
psr = handles.img.coeffs.psr;
soc = handles.img.coeffs.soc;

peaks = xcf * result.ccc + soc * result.autoc + psr * result.psr;
number=str2double(inputdlg('How many particles?'));
if(isnan(result.job.options.modifications.binning))
    result.job.options.modifications.binning = 0;
end;

pickList = tom_os3_returnPicklist(peaks,result.angles ,result.job,number);

plot_pickList(pickList,handles);


combination.xcf = xcf;
combination.psr = psr;
combination.soc = soc;

combinations = handles.filter.combinations;
combinations{length(combinations)+1} = combination;
handles.filter.combinations = combinations;

guidata(hObject, handles);

if(length(combinations) >1)
    set(handles.historySlide,'Enable','on');
    set(handles.historySlide,'Max',length(handles.filter.combinations));
end;
% --- Executes on button press in btnLoad.
function btnLoad_Callback(hObject, eventdata, handles)
% hObject    handle to btnLoad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[FileName,PathName,FilterIndex] = uigetfile;
load([PathName '/' FileName]);

handles.img.coeffs.xcf = combination.xcf;
handles.img.coeffs.psr = combination.psr;
handles.img.coeffs.soc = combination.soc;

set(handles.xcfSlider,'Value',combination.xcf);
set(handles.psrSlider,'Value',combination.psr);
set(handles.socSlider,'Value',combination.soc);

set(handles.xcfEdit,'String',combination.xcf);
set(handles.psrEdit,'String',combination.psr);
set(handles.socEdit,'String',combination.soc);

guidata(hObject, handles);

% --- Executes on button press in btnSave.
function btnSave_Callback(hObject, eventdata, handles)
% hObject    handle to btnSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[FileName,PathName,FilterIndex] = uiputfile;

xcf = handles.img.coeffs.xcf;
psr = handles.img.coeffs.psr;
soc = handles.img.coeffs.soc;

combination.xcf = xcf;
combination.psr = psr;
combination.soc = soc;

save([PathName '/' FileName],'combination');


% --- Executes on button press in rbtn_image.
function rbtn_image_Callback(hObject, eventdata, handles)
% hObject    handle to rbtn_image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rbtn_image

set(handles.rbtn_image,'Value',1);
set(handles.rbtn_peaks,'Value',0);

%%
imageAxis    = handles.imageAxis;
colormap gray;
axes(imageAxis);
imagesc(handles.img.volume');
axis off;

% --- Executes on button press in rbtn_peaks.
function rbtn_peaks_Callback(hObject, eventdata, handles)
% hObject    handle to rbtn_peaks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rbtn_peaks
result = handles.img.result;


set(handles.rbtn_image,'Value',0);
set(handles.rbtn_peaks,'Value',1);

xcf = handles.img.coeffs.xcf;
psr = handles.img.coeffs.psr;
soc = handles.img.coeffs.soc;

peaks = xcf * result.ccc + soc * result.autoc + psr * result.psr; 

imageAxis    = handles.imageAxis;
colormap gray;
axes(imageAxis);
imagesc(peaks');
axis off;


% --- Executes on slider movement.
function historySlide_Callback(hObject, eventdata, handles)
% hObject    handle to historySlide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

if(length(handles.filter.combinations)<2)
    return;
end;

combinations = handles.filter.combinations;
position = get(hObject,'Value')+1;
disp(position);

combinations{ceil(position)}
% --- Executes during object creation, after setting all properties.
function historySlide_CreateFcn(hObject, eventdata, handles)
% hObject    handle to historySlide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in btn_alignRef.
function btn_alignRef_Callback(hObject, eventdata, handles)
% hObject    handle to btn_alignRef (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName PathName] =  uigetfile({'*.em','EM Files'},'Select the aligment reference');

if(~ischar(FileName))
    return;
end;

mask = tom_emread([PathName '/' FileName]);
mask = mask.Value;

options = handles.options;
options.classification.alignRef = mask;
handles.options = options;
guidata(hObject, handles);

% --- Executes on button press in btn_trainingSet.
function btn_trainingSet_Callback(hObject, eventdata, handles)
% hObject    handle to btn_trainingSet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName PathName] =  uigetfile({'*.em','EM Files'},'Select the training set');

if(~ischar(FileName))
    return;
end;

stack = tom_emreadc3([PathName '/' FileName]);
stack = stack.Value;

options = handles.options;
options.classification.trainingStack = stack;
handles.options = options;
guidata(hObject, handles);

% --- Executes on button press in btn_alignMask.
function btn_alignMask_Callback(hObject, eventdata, handles)
% hObject    handle to btn_alignMask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName PathName] = uigetfile({'*.em','EM Files'},'Select mask file');

if(~ischar(FileName))
    return;
end;

mask = tom_emread([PathName '/' FileName]);
mask = mask.Value;

options = handles.options;
options.classification.alignMask = mask;
handles.options = options;
guidata(hObject, handles);



function txt_sigStart_Callback(hObject, eventdata, handles)
% hObject    handle to txt_sigStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_sigStart as text
%        str2double(get(hObject,'String')) returns contents of txt_sigStart as a double
string = get(hObject,'String');

if(strcmp(string,''))
    return;
end;

options = handles.options;
options.classification.sigmaScale(1)= str2double(string);
handles.options = options;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function txt_sigStart_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_sigStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_sigEnd_Callback(hObject, eventdata, handles)
% hObject    handle to txt_sigEnd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_sigEnd as text
%        str2double(get(hObject,'String')) returns contents of txt_sigEnd as a double
string = get(hObject,'String');

if(strcmp(string,''))
    return;
end;

options = handles.options;
options.classification.sigmaScale(2)= str2double(string);
handles.options = options;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function txt_sigEnd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_sigEnd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_sigStep_Callback(hObject, eventdata, handles)
% hObject    handle to txt_sigStep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_sigStep as text
%        str2double(get(hObject,'String')) returns contents of txt_sigStep as a double
string = get(hObject,'String');

if(strcmp(string,''))
    return;
end;

options = handles.options;
options.classification.sigmaScale(3)= str2double(string);
handles.options = options;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function txt_sigStep_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_sigStep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_nrClusters_Callback(hObject, eventdata, handles)
% hObject    handle to txt_nrClusters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_nrClusters as text
%        str2double(get(hObject,'String')) returns contents of txt_nrClusters as a double
string = get(hObject,'String');

if(strcmp(string,''))
    return;
end;

options = handles.options;
options.classification.numberOfClusters= str2double(string);
handles.options = options;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function txt_nrClusters_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_nrClusters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_eigFirst_Callback(hObject, eventdata, handles)
% hObject    handle to txt_eigFirst (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_eigFirst as text
%        str2double(get(hObject,'String')) returns contents of txt_eigFirst as a double
string = get(hObject,'String');

if(strcmp(string,''))
    return;
end;

options = handles.options;
options.classification.eigenStart= str2double(string);
handles.options = options;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function txt_eigFirst_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_eigFirst (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function txt_eigLast_Callback(hObject, eventdata, handles)
% hObject    handle to txt_eigLast (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of txt_eigLast as text
%        str2double(get(hObject,'String')) returns contents of txt_eigLast as a double
string = get(hObject,'String');

if(strcmp(string,''))
    return;
end;

options = handles.options;
options.classification.eigenEnd= str2double(string);
handles.options = options;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function txt_eigLast_CreateFcn(hObject, eventdata, handles)
% hObject    handle to txt_eigLast (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_classifiy.
function btn_classifiy_Callback(hObject, eventdata, handles)
% hObject    handle to btn_classifiy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    

    options = handles.options;
    result = handles.img.result;
    

    
    trainingStack = options.classification.trainingStack;

    if(~isfield(options.classification,'alignMask'))
        options.classification.alignMask = tom_os3_sphereMask(options.classification.alignRef);
    end;
    
    al = questdlg('Is trainingstack aligned?');
    if(strcmp(al,'No'))
        trainingStack = tom_os3_alignStack(trainingStack,options.classification.alignRef,options.classification.alignMask);
    end;
    options.classification.trainingStack = trainingStack;

    xcf = handles.img.coeffs.xcf;
    psr = handles.img.coeffs.psr;
    soc = handles.img.coeffs.soc;

    %linear combination
    peaks = xcf * result.ccc + soc * result.autoc + psr * result.psr;
    number=str2double(INPUTDLG('How many particles?'));
   
    pickList = tom_os3_returnPicklist(peaks,result.angles,result.job,number);
    plot_pickList(pickList,handles);    
    al = questdlg('Binning of?');
    if(strcmp(al,'Yes'))
        result.job.options.modifications.binning =0;
    end;    
    particleStack  = tom_os3_generateParticleStack(pickList,result.job,'');
    particleStack = tom_os3_alignStack(particleStack,options.classification.alignRef,options.classification.alignMask);
    
    maskOn= questdlg('Mask particles?');
    if(strcmp(maskOn,'Yes'))
        
        [FileName PathName] = uigetfile({'*.em','EM Files'},'Select mask file');
        
        m = tom_emreadc3([PathName '/' FileName]);
        m = m.Value;
        
        for i=1:size(options.classification.trainingStack,3)
            options.classification.trainingStack(:,:,i) = options.classification.trainingStack(:,:,i) .* m;
        end;
        
        for i=1:size(particleStack,3)
            particleStack(:,:,i) = particleStack(:,:,i) .* m;
        end;
        
    end;
    plotON = questdlg('Show classification process?');
    plotON = strcmp('Yes',plotON);
    goodFlags = tom_os3_iterativeClassifier(particleStack,options.classification.trainingStack,options,plotON);
    
    pickList = tom_os3_reducePicklist(pickList,goodFlags);
%    align2d = tom_os3_pickList2Align2d(pickList);
    plot_pickList(pickList,handles);
    
%%    
function plot_pickList(pickList,handles)

result = handles.img.result;
    
%select color
button = questdlg('Choose marker color?');
if(strcmp(button,'Yes'))
    colors = uisetcolor;
else
    colors = [0 1 0];
end;

mark = inputdlg('Choose symbol: circle (o : default), square (s), diamond (d)');
if(length(mark) < 1)
    mark = 'o';
else
    mark = mark{1};
end;

axes(handles.imageAxis);
hold on;

for i=1:length(pickList)
    
    pick = pickList{i};
    x = pick.coordinates(1) /2^result.job.options.modifications.binning;
    y = pick.coordinates(2) /2^result.job.options.modifications.binning;
    
    plot(x,y,mark,'Color',colors,'MarkerSize',result.job.templateSize(1)/2^(result.job.options.modifications.binning+1),'LineWidth',3);
    %text(x+5,y+5,num2str(i),'Color',[0 1 0]);
    
end;

hold off;
axis off;


% --------------------------------------------------------------------
function mnu_optimisation_Callback(hObject, eventdata, handles)
% hObject    handle to mnu_optimisation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function mnu_optFilter_Callback(hObject, eventdata, handles)
% hObject    handle to mnu_optFilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


picks = handles.img.picks;
result = handles.img.result;


align2d = [];


for i=1:length(picks)
    
    p = picks{i};
    
    al.filename = result.job.volumeFile;
    al.position.x = p(1);
    al.position.y = p(2);

    al.radius = result.job.templateSize(1)/2^(result.job.options.modifications.binning+1);
    
    
    if(isempty(align2d))
        align2d = al;
    else
        align2d(length(align2d)+1) = al;
    end;
end;


if(~isfield(result.job.options,'filter'))
    result.job.options.filter.numberParticles =str2double(INPUTDLG('How many particles?'));
end;
    
optimum = tom_os3_simulatedAnnealing([0,0,0],[1,1,1],0.1,[0.5,0.5,0.5],10,result,align2d);


