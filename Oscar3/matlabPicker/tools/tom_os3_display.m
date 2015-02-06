function varargout = tom_os3_display(varargin)
% TOM_OS3_DISPLAY M-file for tom_os3_display.fig
%      TOM_OS3_DISPLAY, by itself, creates a new TOM_OS3_DISPLAY or raises the existing
%      singleton*.
%
%      H = TOM_OS3_DISPLAY returns the handle to a new TOM_OS3_DISPLAY or the handle to
%      the existing singleton*.
%
%      TOM_OS3_DISPLAY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TOM_OS3_DISPLAY.M with the given input arguments.
%
%      TOM_OS3_DISPLAY('Property','Value',...) creates a new TOM_OS3_DISPLAY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before tom_os3_display_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to tom_os3_display_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help tom_os3_display

% Last Modified by GUIDE v2.5 30-Nov-2007 15:17:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tom_os3_display_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_os3_display_OutputFcn, ...
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


% --- Executes just before tom_os3_display is made visible.
function tom_os3_display_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to tom_os3_display (see VARARGIN)

% Choose default command line output for tom_os3_display
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes tom_os3_display wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = tom_os3_display_OutputFcn(hObject, eventdata, handles) 
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
plot(x,y,'x','Color',[1 0 0]);
text(x+5,y+5,num2str(handles.img.counter),'Color',[1 0.1 0.1]);
axis off;

% disp([ num2str(handles.img.counter) ' - CCC: ' num2str(handles.img.result.ccc(x,y)) ' FING: ' num2str(handles.img.result.angleCorr(x,y)) ' Prod:' num2str(handles.img.result.autoc(x,y) * handles.img.result.angleCorr(x,y)) ]);
handles.img.counter = handles.img.counter +1;

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
volume = tom_emread(result.job.volumeFile);
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
set(handles.projectionAxis,'Visible','Off');



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
function showCCCOnly_Callback(hObject, eventdata, handles)
% hObject    handle to showCCCOnly (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
result = handles.img.result;
peaks = result.ccc;
number=str2double(INPUTDLG('How many?'));
list = tom_os3_returnPicklist(peaks,result.angles ,result.job,number);


for i=1:length(list)
    
    pick = list{i};
    x = pick.coordinates(1) /2^result.job.options.modifications.binning;
    y = pick.coordinates(2) /2^result.job.options.modifications.binning;
    
    axes(handles.imageAxis);
    hold on;
    plot(x,y,'s','Color',[1 0 0],'MarkerSize',result.job.templateSize(1)/2^(result.job.options.modifications.binning+1)+15,'LineWidth',2);
   % text(x+5,y+5,num2str(i),'Color',[1 0 0]);
    hold off;
    axis off;
    
    axes(handles.correlationAxis);
    hold on;
    plot(x,y,'s','Color',[1 0 0],'MarkerSize',result.job.templateSize(1)/2^(result.job.options.modifications.binning+1),'LineWidth',2);
    %text(x+5,y+5,num2str(i),'Color',[1 0 0]);
    hold off;
    axis off;
end;

% --------------------------------------------------------------------
function showShapeCorrelation_Callback(hObject, eventdata, handles)
% hObject    handle to showShapeCorrelation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
result = handles.img.result;
peaks = handles.img.peaks;
number=str2double(INPUTDLG('How many?'));
if(isnan(result.job.options.modifications.binning))
    result.job.options.modifications.binning = 0;
end;
list = tom_os3_returnPicklist(peaks,result.angles ,result.job,number);


for i=1:length(list)
    
    pick = list{i};
    x = pick.coordinates(1) /2^result.job.options.modifications.binning;
    y = pick.coordinates(2) /2^result.job.options.modifications.binning;
    
    axes(handles.imageAxis);
    hold on;
    plot(x,y,'o','Color',[0 1 0],'MarkerSize',result.job.templateSize(1)/2^(result.job.options.modifications.binning+1)+10,'LineWidth',3);
    %text(x+5,y+5,num2str(i),'Color',[0 1 0]);
    hold off;
    axis off;
    
    axes(handles.correlationAxis);
    hold on;
    plot(x,y,'o','Color',[0 1 0],'MarkerSize',result.job.templateSize(1)/2^(result.job.options.modifications.binning+1),'LineWidth',2);
   % text(x+5,y+5,num2str(i),'Color',[0 1 0]);
    hold off;
    axis off;
end;


% --------------------------------------------------------------------
function adjustImageContrast_Callback(hObject, eventdata, handles)
% hObject    handle to adjustImageContrast (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
imageAxis=handles.imageAxis;
imcontrast(imageAxis);
% --------------------------------------------------------------------
function adjustCorrelationContrast_Callback(hObject, eventdata, handles)
% hObject    handle to adjustCorrelationContrast (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
correlationAxis=handles.correlationAxis;
imcontrast(correlationAxis);

% --------------------------------------------------------------------
function Untitled_6_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Save_CCC_Hits_Callback(hObject, eventdata, handles)
% hObject    handle to Save_CCC_Hits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
number=str2double(INPUTDLG('How many?'));
result = handles.img.result;

list = tom_os3_returnPicklist(result.ccc,result.angles ,result.job.templateSize,number,result.job.options.modifications.binning);

volume = handles.img.volume';
ps = zeros(64,64,number);

for i =1:length(list)
    
    pick = list{i};
    x = pick.coordinates(1) /2^result.job.options.modifications.binning;
    y = pick.coordinates(2) /2^result.job.options.modifications.binning;
    coord= pick.coordinates;
    part = tom_cut_out(handles.img.volume,[x-result.job.templateSize(1)/4 y-result.job.templateSize(2)/4],[result.job.templateSize(1)/2 result.job.templateSize(2)/2],'no-fill');
    part = tom_rotate(tom_norm(part.*mask,'mean0+1std'),-result.job.angleList(pick.angle));
    ps(:,:,i) = part;
    
end;
tom_emwrite([ './ccc' num2str(i) 'Pick.em' ],ps);
align2d=tom_av2_create_alignfromstack('ccc50Pick.em');
try
    align2d(1:length(align2d)).scores  =0;
catch
end
save('ccc50Pick.mat','align2d');
% --------------------------------------------------------------------
function Save_Peak_Hits_Callback(hObject, eventdata, handles)
% hObject    handle to Save_Peak_Hits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
number=str2double(INPUTDLG('How many?'));
result = handles.result;

list = tom_os3_returnPicklist(handles.img.peaks,result.angles,result.job.templateSize*2,number,result.job.options.modifications.binning);


tom_emwrite([ './peak' num2str(i) 'Pick.em' ],ps);
align2d=tom_av2_create_alignfromstack([ './peak' num2str(i) 'Pick.em' ]);
try
    align2d(1:length(align2d)).scores  =0;
catch
end
save([ './peak' num2str(i) 'Pick.mat' ],'align2d');
save(['./hits' num2str(i) '.mat'],'list');
% --------------------------------------------------------------------
function Untitled_9_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_16_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

number=str2double(INPUTDLG('How many?'));
result = handles.img.result;

list = tom_os3_returnPicklist(handles.img.peaks,result.angles ,result.job.templateSize,number,result.job.options.modifications.binning);

volume = handles.volume';
ps = zeros(64,64,number);
mask = (tom_os3_pasteCenter(zeros(64,64),tom_os3_sphereMask(zeros(35,35))));
for i =1:length(list)
    
    pick = list{i};
    x = pick.coordinates(1) /2^result.job.options.modifications.binning;
    y = pick.coordinates(2) /2^result.job.options.modifications.binning;
    coord= pick.coordinates;
    part = tom_cut_out(handles.volume,[x-result.job.templateSize(1)/4 y-result.job.templateSize(2)/4],[result.job.templateSize(1)/2 result.job.templateSize(2)/2],'no-fill');
    part = tom_rotate(tom_norm(part,'mean0+1std'),result.job.angleList(pick.angle));
    ps(:,:,i) = part.*mask;
%     figure;tom_imagesc(part');
    %tom_emwrite([ './peak' num2str(i) 'Pick.em' ],part);
end;


template = tom_emread(result.job.templateFile);
template = tom_bin(template.Value,result.job.options.modifications.binning);


template = tom_os3_pasteCenter(zeros(64,64),template);
template = template .*tom_os3_pasteCenter(zeros(64,64),mask);
figure;
for i=1:number
    d = (template-ps(:,:,i)).^2;
    plot(i,sum(d(:)),'.');
    
end;


% --------------------------------------------------------------------
function loadAngCorr_Callback(hObject, eventdata, handles)
% hObject    handle to loadAngCorr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[AngVolName,AngPathName,FilterIndex] = uigetfile('.em',path);
[ImgName,ImgPathName,FilterIndex] = uigetfile('.em',path);

angVol=tom_emread([AngVolPath '/' AngVolName]);
angVol=angVol.Value;

img=tom_emread([ImgPath '/' ImgPathName]);
img=angVol.Value;


imageAxis    = handles.imageAxis;
colormap gray;
axes(imageAxis);
% --------------------------------------------------------------------
function Untitled_17_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function genAngCorr_Callback(hObject, eventdata, handles)
% hObject    handle to genAngCorr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

job = handles.result.job;

template = tom_emread(job.templateFile);
template = tom_bin(template.Value,job.options.modifications.binning);

%%  generate angular correlation stacks


if(~isfield(handles,'angCorrelation'))
    j.template = template;
    j.angleListOrig = [0:10:359];
    j.options =[];
    j.volume = template;
    [a angTemplateStack] = tom_os3_angularCorrelation(j);
    j.volume = handles.volume;
    [a angVolumeStack] = tom_os3_angularCorrelation(j);

    angCorrelation = tom_os3_corrNEW(angVolumeStack,angTemplateStack);
    do = true;
else
    angCorrelation = handles.angCorrelation;
    do = false;
end;
%%  show projection of stack
correlationAxis = handles.correlationAxis;
axes(correlationAxis);

proj = sum(angCorrelation,3)/size(angCorrelation,3);
imagesc(proj');

if(do)
    handles.img.angCorrelation =  angCorrelation;
    handles.img.angVolumeStack = angVolumeStack;
    handles.img.angTemplateStack = angTemplateStack;
end;
guidata(hObject, handles);
axis off;
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
imagesc(tom_filter(volume',2));
set(get(imageAxis, 'Children'),'ButtonDownFcn','tom_os3_filterSetup(''axisButtonDownCallback'',gcbo,[],guidata(gcbo))')
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
function showAngCorr_Callback(hObject, eventdata, handles)
% hObject    handle to showAngCorr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

result = handles.img.result;
peaks  = handles.img.angCorrelation;
number=str2double(INPUTDLG('How many?'));


radius = result.job.templateSize(1)/2^(result.job.options.modifications.binning+1);

picklist = {};
counter = 1;



while(counter <= number)
    
    [coordinates value] = tom_peakc(peaks);
    
    mask = tom_cylindermask(ones(size(peaks)),radius,0,[coordinates(1) coordinates(2)]);
    mask = mask == 0;
    
    peaks = peaks .*mask;
    
    pick.coordinates = coordinates;
    pick.value       = value;
    pick.templateSize= result.job.templateSize(1);
    picklist{counter} = pick;
    counter = counter +1;
end;


for i=1:length(picklist)
    
    pick = picklist{i};
    x = pick.coordinates(1) ;
    y = pick.coordinates(2) ;
    
    axes(handles.imageAxis);
    hold on;
    plot(x,y,'d','Color',[0 0 1],'MarkerSize',radius,'LineWidth',1);
    text(x+5,y+5,num2str(i),'Color',[0 0 1]);
    hold off;
    axis off;
    
    axes(handles.correlationAxis);
    hold on;
    plot(x,y,'o','Color',[0 0 1],'MarkerSize',radius,'LineWidth',5);
    text(x+5,y+5,num2str(i),'Color',[0 0 1]);
    hold off;
    axis off;
end;



% --------------------------------------------------------------------
function maxAngCorr_Callback(hObject, eventdata, handles)
% hObject    handle to maxAngCorr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
result = handles.img.result;
peaks = zeros(size(handles.angCorrelation,1),size(handles.img.angCorrelation,2),'single');

for x=1:size(peaks,1)
    for y=1:size(peaks,2)

        peaks(x,y) = max(handles.img.angCorrelation(x,y,:));
        
    end;
end;

number=str2double(INPUTDLG('How many?'));

radius = result.job.templateSize(1)/2^(result.job.options.modifications.binning+1);

picklist = {};
counter = 1;

while(counter <= number)
    
    [coordinates value] = tom_peakc(peaks);
    
    mask = tom_cylindermask(ones(size(peaks)),radius,0,[coordinates(1) coordinates(2)]);
    mask = mask == 0;
    
    peaks = peaks .*mask;
    
    pick.coordinates = coordinates;
    pick.value       = value;
    pick.templateSize= result.job.templateSize(1);
    picklist{counter} = pick;
    counter = counter +1;
end;


for i=1:length(picklist)
    
    pick = picklist{i};
    x = pick.coordinates(1) ;
    y = pick.coordinates(2) ;
    
    axes(handles.imageAxis);
    hold on;
    plot(x,y,'d','Color',[0 0 1],'MarkerSize',radius,'LineWidth',1);
    text(x+5,y+5,num2str(i),'Color',[0 0 1]);
    hold off;
    axis off;
    
    axes(handles.correlationAxis);
    hold on;
    plot(x,y,'o','Color',[0 0 1],'MarkerSize',radius,'LineWidth',5);
    text(x+5,y+5,num2str(i),'Color',[0 0 1]);
    hold off;
    axis off;
end;


% --------------------------------------------------------------------
function createStack_Callback(hObject, eventdata, handles)
% hObject    handle to createStack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

result = handles.img.result;
% method = (INPUTDLG('Which method?'));

peaks = get_peaks(result);

number=str2double(INPUTDLG('How many picks?'));


list = tom_os3_returnPicklist(peaks,result.angles,result.job.templateSize*2,number,result.job.options.modifications.binning);

result.job.templateSize = [result.job.templateSize(1) result.job.templateSize(2)];
handles.img.particleStack.stack = tom_os3_generateParticleStack(list,result.job,'none');
handles.img.particleStack.picklist = list;

guidata(hObject, handles);

% --------------------------------------------------------------------
function Untitled_18_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function classifyStack_Callback(hObject, eventdata, handles)
% hObject    handle to classifyStack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



    particleStack = handles.img.particleStack.stack;
    
    %normalize the stack under the mask
    mask = tom_os3_sphereMask(zeros([size(particleStack,1),size(particleStack,2)]));
    for i=1:size(particleStack,3)
        p(:,:,i)= tom_bin(tom_norm(particleStack(:,:,i),'phase',mask),1);
    end;
    particleStack = p;
    
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
    
    axes(handles.correlationAxis);
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
number=str2double(INPUTDLG('How many?'));
if(isnan(result.job.options.modifications.binning))
    result.job.options.modifications.binning = 0;
end;
list = tom_os3_returnPicklist(peaks,result.angles ,result.job,number);


%select color
button = questdlg('Choose marker color?');
if(strcmp(button,'Yes'))
    colors = uisetcolor;
else
    colors = [0 1 0];
end;

for i=1:length(list)
    
    pick = list{i};
    x = pick.coordinates(1) /2^result.job.options.modifications.binning;
    y = pick.coordinates(2) /2^result.job.options.modifications.binning;
    
    axes(handles.imageAxis);
    hold on;
    plot(x,y,'o','Color',colors,'MarkerSize',result.job.templateSize(1)/2^(result.job.options.modifications.binning+1)+10,'LineWidth',3);
    %text(x+5,y+5,num2str(i),'Color',[0 1 0]);
    hold off;
    axis off;
    
    axes(handles.correlationAxis);
    hold on;
    plot(x,y,'o','Color',colors,'MarkerSize',result.job.templateSize(1)/2^(result.job.options.modifications.binning+1),'LineWidth',2);
   % text(x+5,y+5,num2str(i),'Color',[0 1 0]);
    hold off;
    axis off;
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




