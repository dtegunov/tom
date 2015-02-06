function varargout = tom_chooser(varargin)

% TOM_CHOOSER interactively "classifys" particles in tomogram
%   MOTIVELIST_OUT = TOM_CHOOSER(MOTL, VOLUMEFILE, ICLASS, RADIUS, ...
%       DSPCUBMODE, IFILT, INDMAX, STARTIND) 
%
%   TOM_CHOOSER is a GUI designed to classify particles subjectively.
%   A motivelist has to be genertade first e.g. by use of cross-correlation
%   (oscar!). For format of motivelist (motl) see bellow. The location of
%   the potential particle within the tomogram is displayed on a x-y-slice
%   and a gallery of the "particle" is diplayed (top-view). The potential
%   particle can be chosen to be a potetial paticle or not ("Take it" or
%   "Leave it"). The particle get a class number specified by ICLASS.
%
%-->NOTE: In this version a bug in the matlab MONTAGE function is bypassed!
%   This version also works under matlab v7.01 and above.
%
%PARAMETERS
%   MOTL       : MOTIVELIST
%   VOLUMEFILE : Filename of Volume (tomogram!)
%   ICLASS     : Index assigned to "class of particles"
%   RADIUS     : RADIUS of particles (for display) 
%   DSPCUBMODE : mode for av3_dspcub (0,1,2)
%   IFILT      : filtering for display
%   INDMAX     : Number of particles to be exmamined
%   STARTIND   : Index of first tentative particle to be examined
%
%Format of MOTL:
%   The following parameters are stored in the matrix MOTIVELIST of dimension (20, NPARTICLES):
%   column 
%      1         : Cross-Correlation Coefficient
%      2         : x-coordinate in full tomogram
%      3         : y-coordinate in full tomogram
%      4         : peak number
%      5         : running index of tilt series (optional)
%      8         : x-coordinate in full tomogram
%      9         : y-coordinate in full tomogram
%      10        : z-coordinate in full tomogram
%      14        : x-shift in subvolume (AFTER rotation of reference)
%      15        : y-shift in subvolume
%      16        : z-shift in subvolume
%      17        : Phi
%      18        : Psi
%      19        : Theta 
%      20        : class no
%
% EXAMPLE
%   phantom = tom_phantom3d([128 128 128],'tripod');
%   motl = zeros(20,3);
%   motl(8:10,1) = [72,73,71];
%   motl(8:10,2) = [40,41,42];
%   motl(8:10,3) = [42,35,43];
%   motl(4,:) = 1:3;
%   motl_out = tom_chooser(motl,'',1,5,1,2,3,1);
%
%SEE ALSO
%   TOM_PICKER, AV3 package (will be delivered as part of TOM), oscar (C-code)
%
%   10/09/02 FF
%   last change 11/05/04 FF

% TOM_CHOOSER M-file for tom_chooser.fig
%      TOM_CHOOSER, by itself, creates a new TOM_CHOOSER or raises the existing
%      singleton*.
%
%      H = TOM_CHOOSER returns the handle to a new TOM_CHOOSER or the handle to
%      the existing singleton*.
%
%      TOM_CHOOSER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TOM_CHOOSER.M with the given input arguments.
%
%      TOM_CHOOSER('Property','Value',...) creates a new TOM_CHOOSER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before tom_chooser_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to tom_chooser_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help tom_chooser

% Last Modified by GUIDE v2.5 23-Nov-2002 20:39:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @tom_chooser_OpeningFcn, ...
                   'gui_OutputFcn',  @tom_chooser_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT



% --- Executes just before tom_chooser is made visible.
function tom_chooser_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to tom_chooser (see VARARGIN)

% Choose default command line output for tom_chooser
handles.output = hObject;
handles.motl=varargin{1};
handles.volumefile=varargin{2};
handles.header = tom_reademheader(handles.volumefile);
handles.header = handles.header.Header; 
handles.iclass=varargin{3};
handles.radius = varargin{4};
handles.dspcubmode = varargin{5};
handles.ifilt = varargin{6};
if (nargin > 6) 
    if (varargin{7} <= size(handles.motl,2))
        handles.indmax = varargin{7};
    else
        handles.indmax = size(handles.motl,2);
    end;
else
    handles.indmax = size(handles.motl,2);
end;
if (nargin > 7)
    handles.index=varargin{8};
else
    handles.index=1;
end;

% Update handles structure
guidata(hObject, handles);
x=handles.motl(8,1);
y=handles.motl(9,1);
z=handles.motl(10,1);
axes(handles.image);
slice = (tom_emreadc(handles.volumefile,'subregion', [1 1 z-1], [handles.header.Size(1)-1 handles.header.Size(2)-1 0] ));
tom_imagesc(double(slice.Value));
uu = [x x x x-5 x+5];vv = [y-5 y+5 y y y];
line(uu,vv,'LineWidth',1,'Color',[1 0 0]);
particle = tom_emreadc(handles.volumefile,'subregion', [x-handles.radius , y-handles.radius, z-handles.radius], ...
            [2*handles.radius 2*handles.radius 2*handles.radius]);
particle = double(particle.Value);        
particle = double(tom_rotate(particle, [-handles.motl(18,handles.index), ...
              -handles.motl(17,handles.index), -handles.motl(19,handles.index)] ));
disp(['particle ' num2str(handles.index) ': z = ' num2str(z) ...
                    ' theta = ' num2str(handles.motl(19,handles.index)) ...
                    ' psi = ' num2str(handles.motl(18,handles.index))]);     
axes(handles.dspcub);
[meanp, maxp, minp, stdp, varp] = tom_dev(particle,'noinfo');
particle = tom_limit((particle -meanp),-2*stdp,2*stdp);
%av3_dspcub(tom_filter(particle,handles.ifilt,'circ'),handles.dspcubmode); 
tom_imagesc(av3_vol2gallery(tom_filter(particle,handles.ifilt,'circ'),ceil(sqrt(size(particle,1))),handles.dspcubmode));
axes(handles.proj);colormap gray;
particle = sum(particle,3);
tom_imagesc(particle);

% UIWAIT makes tom_chooser wait for user response (see UIRESUME)dspcub
uiwait(gcf);
% uiwait(handles.figure1);



% --- Outputs from this function are returned to the command line.
function varargout = tom_chooser_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure

%tom_emwrite('temp.em',handles.motl);
varargout{1} = handles.motl;
close(gcf);


% --- Executes on button press in take.
function take_Callback(hObject, eventdata, handles)
% hObject    handle to take (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
                      
        handles.motl(20,handles.index) = handles.iclass;
        handles.index=handles.index+1;
        guidata(hObject,handles);
        if handles.index<=handles.indmax
            x=handles.motl(8,handles.index);
            y=handles.motl(9,handles.index);
            z=handles.motl(10,handles.index);
            axes(handles.image); 
            slice = (tom_emreadc(handles.volumefile,'subregion', [1 1 z-1], [handles.header.Size(1)-1 handles.header.Size(2)-1 0] ));
            axes(handles.image);
            tom_imagesc(double(slice.Value));
            uu = [x x x x-5 x+5];vv = [y-5 y+5 y y y];
            line(uu,vv,'LineWidth',1,'Color',[1 0 0]);
            particle = tom_emreadc(handles.volumefile,'subregion', [x-handles.radius , y-handles.radius, z-handles.radius], ...
                        [2*handles.radius 2*handles.radius 2*handles.radius]);
            particle = double(particle.Value);
            particle = double(tom_rotate(particle, [-handles.motl(18,handles.index), ...
                -handles.motl(17,handles.index), -handles.motl(19,handles.index)] ));
            disp(['particle ' num2str(handles.index) ': z = ' num2str(z) ...
                    ' theta = ' num2str(handles.motl(19,handles.index)) ...
                    ' psi = ' num2str(handles.motl(18,handles.index))]);
            axes(handles.dspcub); 
            [meanp, maxp, minp, stdp, varp] = tom_dev(particle,'noinfo');
            particle = tom_limit((particle -meanp),-2*stdp,2*stdp);
            %av3_dspcub(tom_filter(particle,handles.ifilt,'circ'),handles.dspcubmode);
            tom_imagesc(av3_vol2gallery(tom_filter(particle,handles.ifilt,'circ'),ceil(sqrt(size(particle,1))),handles.dspcubmode));
            particle = sum(particle,3);
            axes(handles.proj);
            imagesc(particle');
            guidata(hObject,handles);
        else      
            uiresume(gcf);
        end

% --- Executes on button press in leave.
function leave_Callback(hObject, eventdata, handles)
% hObject    handle to leave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
        
        %handles.motl(20,handles.index) = 0;
        handles.index=handles.index+1;
        guidata(hObject,handles);
        if handles.index<=handles.indmax
            x=handles.motl(8,handles.index);
            y=handles.motl(9,handles.index);
            z=handles.motl(10,handles.index);
            axes(handles.image);
            slice = (tom_emreadc(handles.volumefile,'subregion', [1 1 z-1], [handles.header.Size(1)-1 handles.header.Size(2)-1 0] ));
            tom_imagesc(double(slice.Value));
            uu = [x x x x-5 x+5];vv = [y-5 y+5 y y y];
            line(uu,vv,'LineWidth',1,'Color',[1 0 0]);
            particle = tom_emreadc(handles.volumefile,'subregion', [x-handles.radius , y-handles.radius, z-handles.radius], ...
                        [2*handles.radius 2*handles.radius 2*handles.radius]);
            particle = double(particle.Value);
            particle = double(tom_rotate(particle, [-handles.motl(18,handles.index), ...
                -handles.motl(17,handles.index), -handles.motl(19,handles.index)] ));
            axes(handles.dspcub);
            [meanp, maxp, minp, stdp, varp] = tom_dev(particle,'noinfo');
            particle = tom_limit((particle -meanp),-2*stdp,2*stdp);
            %av3_dspcub(tom_filter(particle,handles.ifilt,'circ'),handles.dspcubmode);
            tom_imagesc(av3_vol2gallery(tom_filter(particle,handles.ifilt,'circ'),ceil(sqrt(size(particle,1))),handles.dspcubmode));
            disp(['particle ' num2str(handles.index) ': z = ' num2str(z) ...
                    ' theta = ' num2str(handles.motl(19,handles.index)) ...
                    ' psi = ' num2str(handles.motl(18,handles.index))]);
            axes(handles.proj);
            particle = sum(particle,3);imagesc(particle');
            guidata(hObject,handles);
        else
            uiresume(gcf);
        end
